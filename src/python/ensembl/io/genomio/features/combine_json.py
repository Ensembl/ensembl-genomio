# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Combines repeat/ncRNA feature JSON documents from split/chunked inputs and performs coordinate liftover."""

import argparse
from collections.abc import Callable, Iterable, Iterator
import hashlib
import json
import logging
from pathlib import Path
import re
from typing import NotRequired, TypeAlias, TypedDict, TypeVar, cast

from ensembl.io.genomio.schemas.json.validate import schema_validator
from ensembl.io.genomio.utils.agp_utils import AgpEntry, build_component_index, lift_range, parse_agp
from ensembl.io.genomio.utils.chunk_utils import get_paths_from_manifest, validate_regex

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

_CHUNK_RE_STRING = r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$"
_SHA256_RE = re.compile(r"^[a-fA-F0-9]{64}$")
_FEATURE_SCHEMA_NAME = "load_features"

JsonValue: TypeAlias = None | bool | int | float | str | list["JsonValue"] | dict[str, "JsonValue"]


class RepeatConsensus(TypedDict):
    repeat_consensus_key: str  # SHA-256 hex
    repeat_name: str
    repeat_class: str
    repeat_type: str
    repeat_consensus: NotRequired[str]


class RepeatFeature(TypedDict):
    seq_region: str
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: str
    repeat_start: int
    repeat_end: int
    repeat_consensus: NotRequired[str | None]
    score: NotRequired[float | None]
    attributes: NotRequired[dict[str, JsonValue]]


class NcrnaFeature(TypedDict):
    seq_region: str
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: str
    biotype: str
    score: NotRequired[float]
    evalue: NotRequired[float]


class _FeatureCoords(TypedDict):
    seq_region: str
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: str


def _normalise_for_comparison(key: str, value: JsonValue) -> JsonValue:
    """
    Returns a normalised copy of a top-level value for equality comparison.

    For 'analysis', ignore run_date so documents produced at different times
    can still be merged as long as the rest of the analysis metadata matches.
    """
    if key == "analysis" and isinstance(value, dict):
        analysis = dict(value)
        analysis.pop("run_date", None)
        return cast(JsonValue, analysis)
    return value


class _TopLevelAccumulator:
    """Accumulates and validates that selected top-level fields are identical across documents."""

    def __init__(self) -> None:
        self._comparison_values: dict[str, JsonValue] = {}
        self._original_values: dict[str, JsonValue] = {}
        self._paths: dict[str, Path] = {}

    def require_same(self, key: str, value: JsonValue, path: Path) -> None:
        """
        Registers a top-level field value and ensures consistency.

        Args:
            key: Top-level JSON field name (e.g. "analysis", "source").
            value: Parsed JSON value associated with the field.
            path: Path of the document providing this value (used for error reporting).

        Raises:
            ValueError:
                If the same key has already been seen with a different value
                in a previous document.
        """
        date_removed_value = _normalise_for_comparison(key, value)

        if key not in self._comparison_values:
            self._comparison_values[key] = date_removed_value
            self._original_values[key] = value
            self._paths[key] = path
            return

        prev = self._comparison_values[key]
        if prev != date_removed_value:
            prev_path = self._paths[key]
            raise ValueError(
                f"Top-level '{key}' differs between inputs:\n"
                f"  - {prev_path}\n"
                f"  - {path}\n"
                "This script currently requires them to be identical."
            )

    def get_required(self, key: str) -> JsonValue:
        """
        Returns the stored value for a required top-level field.

        Args:
            key: Field name that must have been previously registered.

        Returns:
            The JSON value previously recorded for the field.

        Raises:
            ValueError:
                If the key was never observed (e.g. empty manifest or
                required field missing from all inputs).
        """
        if key not in self._original_values:
            raise ValueError(f"Missing required top-level '{key}' (empty manifest?)")
        return self._original_values[key]


Feature = TypeVar("Feature", bound=_FeatureCoords)

CoerceFeatureFn = Callable[[JsonValue, Path], Feature]


def _get_agp_entry_for_range(
    parts: list[AgpEntry], start: int, end: int, component_id: str, path: Path
) -> AgpEntry:
    """
    Retrieves the unique AGP part whose component span fully contains a given range.

    Args:
        parts: AGP entries for a single component_id.
        start: Feature start on the component (1-based, inclusive).
        end: Feature end on the component (1-based, inclusive).
        component_id: Component identifier being lifted (used for error messages).
        path: Path of the source file being processed (used for error messages).

    Returns:
        The single AGP entry whose ``part_start``/``part_end`` contains the range.

    Raises:
        ValueError: If no entry contains the range, or if multiple entries match (ambiguous mapping).
    """
    matches = [pt for pt in parts if pt.part_start <= start and end <= pt.part_end]
    if not matches:
        spans = ", ".join(f"{pt.part_start}-{pt.part_end}@{pt.record}" for pt in parts)
        raise ValueError(
            f"Range {start}-{end} for component '{component_id}' does not fit within any AGP span ({path}). "
            f"Spans: {spans}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous AGP mapping for component '{component_id}' range {start}-{end} ({path}); "
            f"matched {len(matches)} entries."
        )
    return matches[0]


def _iterate_validated_documents(
    json_paths: Iterable[Path], validate: bool = True
) -> Iterator[tuple[Path, dict[str, JsonValue]]]:
    """
    Yields (path, document) pairs after schema-validating and loading each JSON file.

    Args:
        json_paths: Iterable of input JSON paths (optionally gzipped).
        validate: If False, skips schema validation.

    Yields:
        Tuples of (path, parsed_document).

    Raises:
        ValueError:
            If schema validation fails for any input file (when ``validate=True``).
    """
    for p in json_paths:
        if validate:
            schema_validator(json_file=p, json_schema=_FEATURE_SCHEMA_NAME)
        yield p, _load_json_document(p)


def _load_json_document(path: Path) -> dict[str, JsonValue]:
    """Loads a single JSON object from (optionally gzipped) file."""
    with open_gz_file(path) as fh:
        raw = fh.read()

    text = raw.decode("utf-8") if isinstance(raw, bytes) else raw
    parsed: JsonValue = json.loads(text)

    return cast(dict[str, JsonValue], parsed)


def _write_and_validate(out_json: Path, combined_json: dict[str, JsonValue]) -> None:
    """
    Writes a combined JSON document to disk and validates it against the schema.

    Args:
        out_json: Destination file path for the combined JSON document.
        combined_json: Fully constructed JSON object to serialise.

    Raises:
        ValueError:
            If schema validation fails.
    """
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(
        json.dumps(combined_json, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    schema_validator(out_json, json_schema=_FEATURE_SCHEMA_NAME)


def _detect_load_type(document: dict[str, JsonValue], path: Path) -> str:
    """Detects the load type (repeat vs ncRNA) of a JSON document based on presence of top-level keys."""
    if "repeat_features" in document:
        return "repeat"
    if "ncrna_features" in document and "ncrna_tool" in document:
        return "ncrna"
    raise ValueError(
        f"Cannot auto-detect schema kind for {path}. Expected one of top-level keys: "
        "'repeat_features' or 'ncrna_features'."
    )


def _is_sha256_hex(s: str) -> bool:
    return bool(_SHA256_RE.fullmatch(s))


def _norm_consensus(seq: str) -> str:
    return "".join(seq.split()).upper()


def _compute_consensus_sha256(rn: str, rc_class: str, rt: str, rc_seq: str | None) -> str:
    norm = _norm_consensus(rc_seq) if rc_seq else ""
    payload = "\t".join([rn, rc_class, rt, norm])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _coerce_repeat_consensus(obj: JsonValue, source_path: Path) -> RepeatConsensus:
    """Validates and casts an arbitrary object to RepeatConsensus."""
    rc = cast(RepeatConsensus, obj)
    expected = _compute_consensus_sha256(
        rc["repeat_name"],
        rc["repeat_class"],
        rc["repeat_type"],
        rc.get("repeat_consensus"),
    )
    if rc["repeat_consensus_key"].lower() != expected:
        raise ValueError(
            f"repeat_consensus_key mismatch for {rc['repeat_name']}|{rc['repeat_class']}|{rc['repeat_type']} "
            f"({source_path}). expected={expected} got={rc['repeat_consensus_key']}"
        )
    return rc


def _coerce_repeat_feature(obj: JsonValue, source_path: Path) -> RepeatFeature:
    """Validates and casts an arbitrary object to RepeatFeature (minimal fields)."""
    feat = cast(RepeatFeature, obj)
    if feat["seq_region_start"] > feat["seq_region_end"]:
        raise ValueError(
            f"seq_region_start > seq_region_end for {feat['seq_region']}: "
            f"{feat['seq_region_start']}>{feat['seq_region_end']} ({source_path})"
        )
    return feat


def _coerce_ncrna_feature(obj: JsonValue, source_path: Path) -> NcrnaFeature:
    """Validates and casts an arbitrary object to NcrnaFeature."""
    feat = cast(NcrnaFeature, obj)
    if feat["seq_region_start"] > feat["seq_region_end"]:
        raise ValueError(
            f"seq_region_start > seq_region_end for {feat['seq_region']}: "
            f"{feat['seq_region_start']}>{feat['seq_region_end']} ({source_path})"
        )
    return feat


def _merge_repeat_consensus(
    combined: dict[str, RepeatConsensus],
    incoming: Iterable[JsonValue],
    source_path: Path,
) -> None:
    """Merges repeat_consensus entries."""
    for raw in incoming:
        rc = _coerce_repeat_consensus(raw, source_path)
        key = rc["repeat_consensus_key"].lower()
        combined[key] = rc


def _feature_consensus_key(feature: RepeatFeature, source_path: Path) -> str | None:
    """Validates and returns consensus key from a feature."""
    ref = feature.get("repeat_consensus")
    if ref is None:
        return None
    if not isinstance(ref, str) or not _is_sha256_hex(ref):
        raise ValueError(f"repeat_feature.repeat_consensus must be a 64-char SHA-256 hex ({source_path})")
    return ref.lower()


def _lift_feature_coords(
    feature: Feature,
    chunk_re: re.Pattern[str],
    agp_by_component: dict[str, list[AgpEntry]] | None,
    allow_revcomp: bool,
    source_path: Path,
) -> Feature:
    """
    Constructs a modified copy of a feature with lifted seq_region and coordinates.

    Uses AGP-driven lifting when ``agp_by_component`` is provided; otherwise uses seq_region-driven
    lifting via ``chunk_re``. If the feature does not appear chunked in seq_region-driven mode, it
    is returned unchanged.

    Args:
        feature: Feature mapping containing seq_region and coordinates to lift.
        chunk_re: Regex identifying chunked seq_region IDs; must include named groups ``base`` and ``start``.
        agp_by_component: Optional mapping from AGP component_id to its AGP entries. If provided, AGP-driven
            lifting is used.
        allow_revcomp: Whether to allow reverse-oriented AGP entries (orientation '-') when lifting.
        source_path: Path of the source JSON file (used for error messages).

    Returns:
        A copy of the input feature with updated seq_region / start / end / strand as needed,
        or the original feature if no lifting is required.

    Raises:
        ValueError: If seq_region_start > seq_region_end, if seq_region_strand is invalid, or if AGP mapping
            is ambiguous/invalid.
        KeyError: If AGP-driven lifting is requested but feature.seq_region is not present in the AGP index.
    """
    name = feature["seq_region"]
    sr_start = feature["seq_region_start"]
    sr_end = feature["seq_region_end"]
    strand = feature["seq_region_strand"]

    if sr_start > sr_end:
        raise ValueError(f"seq_region_start > seq_region_end for {name}: {sr_start}>{sr_end} ({source_path})")

    out: dict[str, JsonValue] = cast(dict[str, JsonValue], dict(feature))

    if agp_by_component is not None:
        parts = agp_by_component.get(name)
        if not parts:
            raise KeyError(f"seq_region '{name}' not found as an AGP component_id ({source_path})")

        part = _get_agp_entry_for_range(parts, sr_start, sr_end, component_id=name, path=source_path)
        obj_id, new_start, new_end = lift_range(part, sr_start, sr_end, allow_revcomp)

        out["seq_region"] = obj_id
        out["seq_region_start"] = new_start
        out["seq_region_end"] = new_end

        if part.orientation == "-":
            if strand == "+":
                out["seq_region_strand"] = "-"
            elif strand == "-":
                out["seq_region_strand"] = "+"
            else:
                raise ValueError(f"Invalid strand value: {strand}")

        return cast(Feature, out)

    m = chunk_re.match(name)
    if not m:
        return feature  # unchunked; nothing to do

    base = m.group("base")
    offset = int(m.group("start")) - 1

    out["seq_region"] = base
    out["seq_region_start"] = sr_start + offset
    out["seq_region_end"] = sr_end + offset
    return cast(Feature, out)


def _combine_feature_docs(
    documents: Iterable[tuple[Path, dict[str, JsonValue]]],
    feature_list_key: str,
    coerce_feature: CoerceFeatureFn,
    chunk_re: re.Pattern[str],
    agp_by_component: dict[str, list[AgpEntry]] | None,
    allow_revcomp: bool,
    required_top_level_keys: list[str],
) -> tuple[dict[str, JsonValue], list[Feature], int]:
    """
    Generic combination engine for feature-based JSON documents.

    This function implements the shared logic for combining multiple
    schema-valid JSON documents that contain:
        - One or more required top-level metadata fields that must be identical across all inputs
        (e.g. ``analysis``, ``source``, ``ncrna_tool``).
        - A list of feature objects under a known key (e.g. ``repeat_features`` or ``ncrna_features``).

    Args:
        documents: Iterable of ``(path, document)`` pairs.
        feature_list_key: Top-level key containing the list of feature objects.
        coerce_feature: Function that validates and casts raw feature JSON
            objects into the expected TypedDict form.
        chunk_re: Regex identifying chunked seq_region identifiers when
            AGP-driven lifting is not used.
        agp_by_component: Optional AGP component index enabling AGP-driven lifting.
        allow_revcomp: Whether to allow reverse-orientation AGP lifting.
        required_top_level_keys: List of top-level keys that must be identical
            across all input documents.

    Returns:
        A tuple containing:
            - A dict of required top-level metadata values.
            - The combined list of lifted features.
            - The number of features read across all inputs.

    Raises:
        ValueError:
            If:
                - required top-level fields differ between inputs,
                - no documents are provided,
                - a feature fails coercion or has invalid coordinates,
                - header-driven lifting fails due to malformed chunk identifiers.
        KeyError:
            If AGP-driven lifting is requested but a feature seq_region
            is not present in the AGP index.
    """
    acc = _TopLevelAccumulator()
    combined_features: list[Feature] = []
    nr_features = 0

    any_docs = False
    for p, document in documents:
        any_docs = True

        for k in required_top_level_keys:
            acc.require_same(key=k, value=document[k], path=p)

        raw_features = cast(list[JsonValue], document[feature_list_key])
        for raw in raw_features:
            nr_features += 1
            feat = coerce_feature(raw, p)
            lifted = _lift_feature_coords(
                feature=feat,
                chunk_re=chunk_re,
                agp_by_component=agp_by_component,
                allow_revcomp=allow_revcomp,
                source_path=p,
            )
            combined_features.append(lifted)

    if not any_docs:
        raise ValueError("No JSON files were read from the manifest (empty manifest?)")

    top_level: dict[str, JsonValue] = {k: acc.get_required(k) for k in required_top_level_keys}
    return top_level, combined_features, nr_features


def _combine_repeat_json_paths(
    json_paths: list[Path],
    out_json: Path,
    chunk_re: re.Pattern[str],
    agp_by_component: dict[str, list[AgpEntry]] | None,
    allow_revcomp: bool,
) -> None:
    """
    Combines JSON documents containing repeat features and consensus, with coordinate liftover.

    Args:
        json_paths: Input JSON file paths (optionally gzipped).
        out_json: Output path for the combined JSON document.
        chunk_re: Regex identifying chunked seq_region IDs for seq_region-driven lifting.
        agp_by_component: Optional AGP component index enabling AGP-driven lifting.
        allow_revcomp: Whether to allow reverse-oriented AGP entries when lifting.

    Raises:
        ValueError:
            If:
                - input files differ in required top-level fields (``analysis``, ``source``, ``ncrna_tool``),
                - no input files are provided,
                - a feature has invalid coordinates (e.g. start > end),
                - coordinate lifting fails due to invalid or ambiguous mapping.
        KeyError:
            If AGP-driven lifting is requested but a feature.seq_region is not present in the AGP
            component index.
    """
    consensus_by_key: dict[str, RepeatConsensus] = {}

    for p, document in _iterate_validated_documents(json_paths):
        incoming = cast(list[JsonValue], document.get("repeat_consensus", []))
        _merge_repeat_consensus(consensus_by_key, incoming, p)

    top_level, features, nr_features = _combine_feature_docs(
        documents=_iterate_validated_documents(json_paths, validate=False),
        feature_list_key="repeat_features",
        coerce_feature=_coerce_repeat_feature,
        chunk_re=chunk_re,
        agp_by_component=agp_by_component,
        allow_revcomp=allow_revcomp,
        required_top_level_keys=["analysis", "source"],
    )

    repeat_features = cast(list[RepeatFeature], features)

    missing: set[str] = set()
    for feat in repeat_features:
        key = _feature_consensus_key(feat, source_path=out_json)
        if key is not None and key not in consensus_by_key:
            missing.add(key)

    if missing:
        missing_list = ", ".join(sorted(missing)[:20])
        extra = "" if len(missing) <= 20 else f" (+{len(missing) - 20} more)"
        raise ValueError(
            "repeat_features reference repeat_consensus key(s) not present in repeat_consensus: "
            f"{missing_list}{extra}"
        )

    combined_json: dict[str, JsonValue] = {
        "analysis": top_level["analysis"],
        "source": top_level["source"],
        "repeat_consensus": cast(JsonValue, list(consensus_by_key.values())),
        "repeat_features": cast(JsonValue, repeat_features),
    }

    _write_and_validate(out_json, combined_json)
    logging.info(f"Read {nr_features} repeat_feature(s); wrote {len(features)} repeat_feature(s)")


def _combine_ncrna_json_paths(
    json_paths: list[Path],
    out_json: Path,
    chunk_re: re.Pattern[str],
    agp_by_component: dict[str, list[AgpEntry]] | None,
    allow_revcomp: bool,
) -> None:
    """
    Combines JSON documents containing ncRNA features for a single tool, with coordinate liftover.

    Args:
        json_paths: Input JSON file paths (optionally gzipped). All files must conform to the
            ``load_features`` schema and contain ``ncrna_features`` and ``ncrna_tool`` top-level keys.
        out_json: Output path for the combined JSON document.
        chunk_re: Regex identifying chunked seq_region IDs for seq_region-driven coordinate lifting.
            Must define named groups ``base`` and ``start``.
        agp_by_component: Optional AGP component index enabling AGP-driven lifting. If provided,
            feature.seq_region is treated as an AGP component_id and coordinates are lifted into AGP
            object coordinates.
        allow_revcomp: Whether to allow reverse-oriented AGP entries when lifting.
            If True and an AGP entry has orientation '-', feature strand will be flipped.

    Raises:
        ValueError:
            If:
                - input files differ in required top-level fields,
                - no input files are provided,
                - a feature has invalid coordinates (e.g. start > end),
                - coordinate lifting fails due to ambiguous or invalid AGP spans.
        KeyError:
            If AGP-driven lifting is requested but a feature.seq_region
            is not present in the AGP component index.
    """
    top_level, features, nr_features = _combine_feature_docs(
        documents=_iterate_validated_documents(json_paths),
        feature_list_key="ncrna_features",
        coerce_feature=_coerce_ncrna_feature,
        chunk_re=chunk_re,
        agp_by_component=agp_by_component,
        allow_revcomp=allow_revcomp,
        required_top_level_keys=["analysis", "source", "ncrna_tool"],
    )

    combined_json: dict[str, JsonValue] = {
        "analysis": top_level["analysis"],
        "source": top_level["source"],
        "ncrna_tool": top_level["ncrna_tool"],
        "ncrna_features": cast(JsonValue, cast(list[NcrnaFeature], features)),
    }

    _write_and_validate(out_json, combined_json)
    logging.info(f"Read {nr_features} ncrna_feature(s); wrote {len(features)} ncrna_feature(s)")


def combine_feature_json(
    json_manifest: Path,
    out_json: Path,
    chunk_re: re.Pattern[str] = re.compile(_CHUNK_RE_STRING),
    agp_file: Path | None = None,
    allow_revcomp: bool = False,
) -> None:
    """
    Combines JSON documents from split/chunked inputs and performs coordinate liftover.

    Supports load_features schema documents containing either repeat features with (optional) consensus,
    or ncRNA features for a single tool.

    Liftover modes:
    1) AGP-driven (recommended):
        - Treats feature.seq_region as AGP component_id (part_id).
        - Lifts seq_region_start/end into AGP object coordinates (record_start/end).
        - Replaces seq_region with AGP object id (record).
        - If AGP orientation is '-', flips seq_region_strand (requires --allow-revcomp).

    2) Header-driven (no AGP):
        - Parses seq_region with ``chunk-id-regex``
        - Replaces seq_region with base id and shift coordinates by (chunk_start - 1).

    Merging:
    - Features are concatenated after lifting.
    - Top-level metadata must be identical across input files.

    Args:
        json_manifest: Path to a manifest file containing one JSON file path
            per line. Files may be gzipped.
        out_json: Destination path for the combined JSON output.
        chunk_re: Compiled regex used to identify chunked seq_region identifiers
            and extract coordinate offsets when AGP-driven lifting is not used.
            Must define named groups ``base`` and ``start``.
        agp_file: Optional AGP file enabling AGP-driven coordinate lifting.
            If provided, feature.seq_region values are interpreted as AGP
            component IDs.
        allow_revcomp: Whether to allow reverse-oriented AGP entries when
            lifting coordinates. If True and an AGP entry has orientation '-',
            feature strand will be flipped.

    Raises:
        ValueError:
            If:
                - the manifest is empty,
                - mixed load types (repeat vs ncRNA) are detected in the manifest,
                - required top-level metadata fields differ between inputs,
                - consensus validation fails (repeat mode),
                - features reference missing repeat_consensus keys (repeat mode),
                - feature coordinates are invalid or coordinate lifting is ambiguous.
        KeyError:
            If AGP-driven lifting is requested but a feature.seq_region is not present in the AGP
            component index.
    """
    json_paths = get_paths_from_manifest(json_manifest)
    if not json_paths:
        raise ValueError("No JSON files were read from the manifest (empty manifest?)")

    agp_by_component: dict[str, list[AgpEntry]] | None = None
    if agp_file is not None:
        agp_records = parse_agp(agp_file, allow_revcomp)
        agp_by_component = build_component_index(agp_records)

    # Detect from first document, then enforce all documents match.
    first_doc = _load_json_document(json_paths[0])
    load_type = _detect_load_type(first_doc, json_paths[0])

    mismatched: list[Path] = []
    for p in json_paths[1:]:
        doc = _load_json_document(p)
        if _detect_load_type(doc, p) != load_type:
            mismatched.append(p)

    if mismatched:
        paths = "\n".join(f"  - {p}" for p in mismatched[:20])
        extra = "" if len(mismatched) <= 20 else f"\n  (+{len(mismatched) - 20} more)"
        raise ValueError(
            f"Mixed load types detected in manifest {json_manifest}. "
            f"First file is '{load_type}', but these differ:\n{paths}{extra}"
        )

    if load_type == "repeat":
        _combine_repeat_json_paths(
            json_paths,
            out_json,
            chunk_re=chunk_re,
            agp_by_component=agp_by_component,
            allow_revcomp=allow_revcomp,
        )
    else:
        _combine_ncrna_json_paths(
            json_paths,
            out_json,
            chunk_re=chunk_re,
            agp_by_component=agp_by_component,
            allow_revcomp=allow_revcomp,
        )

    logging.info(f"Wrote combined file to {out_json}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parses command-line arguments for the feature JSON combining CLI.

    Args:
        argv: Optional argument vector (excluding program name). If None, arguments are read from
        ``sys.argv``.

    Returns:
        Parsed argparse namespace with validated options and logging configuration applied.
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--json-manifest",
        metavar="TXT",
        required=True,
        help=(
            "Manifest file containing paths of schema-valid JSON files to be combined (files may "
            "be gzipped)."
        ),
    )
    parser.add_argument_dst_path(
        "--out-json",
        metavar="JSON",
        required=True,
        help="Path for combined JSON output.",
    )
    parser.add_argument_src_path(
        "--agp-file",
        metavar="AGP",
        default=argparse.SUPPRESS,
        help="Optional AGP file; if provided, coordinate lifting uses AGP component mapping.",
    )
    parser.add_argument(
        "--allow-revcomp",
        action="store_true",
        help="Allow reverse-orientation lifting if AGP orientation is '-' (also flips feature strand).",
    )
    parser.add_argument(
        "--chunk-id-regex",
        type=str,
        metavar="REGEX",
        default=_CHUNK_RE_STRING,
        help=(
            "Regex used to identify chunked sequence IDs and extract coordinates when no AGP is provided. "
            "Must define named groups 'base' and 'start'."
        ),
    )
    parser.add_log_arguments()

    args = parser.parse_args(argv)
    init_logging_with_args(args)
    return args


def main(argv: list[str] | None = None) -> None:
    """Entry point for the feature JSON combining CLI."""
    args = parse_args(argv)
    chunk_re = validate_regex(args.chunk_id_regex)
    try:
        combine_feature_json(
            json_manifest=args.json_manifest,
            out_json=args.out_json,
            agp_file=getattr(args, "agp_file", None),
            allow_revcomp=args.allow_revcomp,
            chunk_re=chunk_re,
        )
    except Exception:
        logging.exception(f"Error combining feature JSON from files in {args.json_manifest}")
        raise


if __name__ == "__main__":
    main()
