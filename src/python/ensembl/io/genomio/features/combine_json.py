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
from collections.abc import Iterable
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

_MD5_RE = re.compile(r"^[a-fA-F0-9]{32}$")
_FEATURE_SCHEMA_NAME = "load_features"

JsonValue: TypeAlias = None | bool | int | float | str | list["JsonValue"] | dict[str, "JsonValue"]


class RepeatConsensus(TypedDict):
    repeat_consensus_key: str  # md5hex
    repeat_name: str
    repeat_class: str
    repeat_type: str
    repeat_consensus: NotRequired[str]


class RepeatFeature(TypedDict):
    seq_region: str
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: int
    repeat_start: int
    repeat_end: int
    repeat_consensus: NotRequired[str | None]
    score: NotRequired[float | None]
    attributes: NotRequired[dict[str, JsonValue]]


class NcrnaFeature(TypedDict):
    seq_region: str
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: int
    biotype: str
    score: NotRequired[float]
    evalue: NotRequired[float]


class _FeatureCoords(TypedDict):
    seq_region: str
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: int


Feature = TypeVar("Feature", bound=_FeatureCoords)


def _get_agp_entry_for_range(
    parts: list[AgpEntry], start: int, end: int, *, component_id: str, path: Path
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


def _load_json_document(path: Path) -> dict[str, JsonValue]:
    """Loads a single JSON object from (optionally gzipped) file."""
    with open_gz_file(path) as fh:
        raw = fh.read()

    text = raw.decode("utf-8") if isinstance(raw, bytes) else raw
    parsed: JsonValue = json.loads(text)

    return cast(dict[str, JsonValue], parsed)


def _schema_validate_file(path: Path) -> None:
    """Validates a JSON file against a repository schema using schema_validator."""
    if path.suffix != ".gz":
        schema_validator(json_file=path, json_schema=_FEATURE_SCHEMA_NAME)
        return

    tmp = path.with_suffix("")
    if tmp.exists():
        tmp = path.with_name(f"{tmp.name}.tmp_unzipped_for_schema_validation.json")

    logging.debug(f"Decompressing {path} -> {tmp} for schema validation")

    with open_gz_file(path) as in_fh, tmp.open("w", encoding="utf-8") as out_fh:
        raw = in_fh.read()
        text = raw.decode("utf-8") if isinstance(raw, bytes) else raw
        out_fh.write(text)

    try:
        schema_validator(json_file=tmp, json_schema=_FEATURE_SCHEMA_NAME)
    finally:
        try:
            tmp.unlink(missing_ok=True)
        except Exception:
            logging.warning(f"Failed to remove temp validation file {tmp}", exc_info=True)


def _assert_same_top_level(
    meta_name: str,
    a: dict[str, JsonValue],
    b: dict[str, JsonValue],
    path_a: Path,
    path_b: Path,
) -> None:
    """Ensures top-level metadata objects are identical."""
    if a != b:
        raise ValueError(
            f"Top-level '{meta_name}' differs between inputs:\n"
            f"  - {path_a}\n"
            f"  - {path_b}\n"
            "This script currently requires them to be identical."
        )


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


def _is_md5_hex(s: str) -> bool:
    return bool(_MD5_RE.fullmatch(s))


def _norm_consensus(seq: str) -> str:
    return "".join(seq.split()).upper()


def _compute_consensus_md5(rn: str, rc_class: str, rt: str, rc_seq: str | None) -> str:
    norm = _norm_consensus(rc_seq) if rc_seq else ""
    payload = "\t".join([rn, rc_class, rt, norm])
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def _coerce_repeat_consensus(obj: JsonValue, source_path: Path) -> RepeatConsensus:
    """Validates and casts an arbitrary object to RepeatConsensus."""
    rc = cast(RepeatConsensus, obj)
    expected = _compute_consensus_md5(
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
    """Merges repeat_consensus entries; duplicates must be identical."""
    for raw in incoming:
        rc = _coerce_repeat_consensus(raw, source_path=source_path)
        key = rc["repeat_consensus_key"].lower()

        existing = combined.get(key)
        if existing is None:
            combined[key] = rc
            continue
        if existing != rc:
            raise ValueError(
                f"Conflicting repeat_consensus for key={key} ({source_path}). "
                "Entries differ; cannot safely merge."
            )


def _feature_consensus_key(feature: RepeatFeature, source_path: Path) -> str | None:
    """Validates and returns consensus key from a feature."""
    ref = feature.get("repeat_consensus")
    if ref is None:
        return None
    if not isinstance(ref, str) or not _is_md5_hex(ref):
        raise ValueError(f"repeat_feature.repeat_consensus must be a 32-char hex md5 ({source_path})")
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
        ValueError: If seq_region_start > seq_region_end, or if AGP mapping is ambiguous/invalid.
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
            out["seq_region_strand"] = -strand

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
        ValueError: If:
            - input files differ in top-level analysis or source,
            - consensus entries conflict,
            - features reference missing consensus keys,
            - coordinates are invalid or ambiguous,
            - no input files are provided.
        KeyError: If AGP-driven lifting is requested but a feature seq_region is not present in the AGP index.
    """

    combined_analysis: dict[str, JsonValue] | None = None
    combined_source: dict[str, JsonValue] | None = None
    analysis_path: Path | None = None
    source_path: Path | None = None

    consensus_by_key: dict[str, RepeatConsensus] = {}
    combined_features: list[RepeatFeature] = []

    n_features_in = 0

    for p in json_paths:
        _schema_validate_file(p)
        document = _load_json_document(p)

        analysis = cast(dict[str, JsonValue], document["analysis"])
        source = cast(dict[str, JsonValue], document["source"])

        if combined_analysis is None:
            combined_analysis = analysis
            analysis_path = p
        else:
            _assert_same_top_level(
                "analysis", combined_analysis, analysis, path_a=analysis_path or p, path_b=p
            )

        if combined_source is None:
            combined_source = source
            source_path = p
        else:
            _assert_same_top_level("source", combined_source, source, path_a=source_path or p, path_b=p)

        incoming_consensus = cast(list[JsonValue], document.get("repeat_consensus", []))
        _merge_repeat_consensus(consensus_by_key, incoming_consensus, p)

        features = cast(list[JsonValue], document["repeat_features"])
        for raw_feat in features:
            n_features_in += 1
            feat = _coerce_repeat_feature(raw_feat, source_path=p)
            lifted = _lift_feature_coords(
                feature=feat,
                chunk_re=chunk_re,
                agp_by_component=agp_by_component,
                allow_revcomp=allow_revcomp,
                source_path=p,
            )
            combined_features.append(lifted)

    if combined_analysis is None or combined_source is None:
        raise ValueError("No JSON files were read from the manifest (empty manifest?)")

    missing: set[str] = set()
    for feat in combined_features:
        key = _feature_consensus_key(feat, source_path=out_json)
        if key is None:
            continue
        if key not in consensus_by_key:
            missing.add(key)

    if missing:
        missing_list = ", ".join(sorted(missing)[:20])
        extra = "" if len(missing) <= 20 else f" (+{len(missing) - 20} more)"
        raise ValueError(
            "repeat_features reference repeat_consensus key(s) not present in repeat_consensus: "
            f"{missing_list}{extra}"
        )

    combined_json: dict[str, JsonValue] = {
        "analysis": cast(JsonValue, combined_analysis),
        "source": cast(JsonValue, combined_source),
        "repeat_consensus": cast(JsonValue, list(consensus_by_key.values())),
        "repeat_features": cast(JsonValue, combined_features),
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w", encoding="utf-8") as out_fh:
        json.dump(combined_json, out_fh, ensure_ascii=False, indent=2)
        out_fh.write("\n")

    schema_validator(out_json, json_schema=_FEATURE_SCHEMA_NAME)

    logging.info(f"Read {n_features_in} repeat_feature(s); wrote {len(combined_features)} repeat_feature(s)")


def _combine_ncrna_json_paths(
    json_paths: list[Path],
    out_json: Path,
    *,
    chunk_re: re.Pattern[str],
    agp_by_component: dict[str, list[AgpEntry]] | None,
    allow_revcomp: bool,
) -> None:
    """
    Combines JSON documents containing ncRNA features for a single tool, with coordinate liftover.

    Args:
        json_paths: Input JSON file paths (optionally gzipped).
        out_json: Output path for the combined JSON document.
        chunk_re: Regex identifying chunked seq_region IDs for seq_region-driven lifting.
        agp_by_component: Optional AGP component index enabling AGP-driven lifting.
        allow_revcomp: Whether to allow reverse-oriented AGP entries when lifting.

    Raises:
        ValueError: If top-level metadata differs between inputs, if ncrna_tool differs between inputs,
            if no inputs are provided, or if any feature has invalid coordinates.
        KeyError: If AGP-driven lifting is requested but a feature seq_region is not present in the AGP index.
    """
    combined_analysis: dict[str, JsonValue] | None = None
    combined_source: dict[str, JsonValue] | None = None
    combined_tool: str | None = None

    analysis_path: Path | None = None
    source_path: Path | None = None
    tool_path: Path | None = None

    combined_features: list[NcrnaFeature] = []
    n_features_in = 0

    for p in json_paths:
        _schema_validate_file(p)
        document = _load_json_document(p)

        analysis = cast(dict[str, JsonValue], document["analysis"])
        source = cast(dict[str, JsonValue], document["source"])
        tool = cast(str, document["ncrna_tool"])

        if combined_analysis is None:
            combined_analysis = analysis
            analysis_path = p
        else:
            _assert_same_top_level(
                "analysis", combined_analysis, analysis, path_a=analysis_path or p, path_b=p
            )

        if combined_source is None:
            combined_source = source
            source_path = p
        else:
            _assert_same_top_level("source", combined_source, source, path_a=source_path or p, path_b=p)

        if combined_tool is None:
            combined_tool = tool
            tool_path = p
        else:
            if combined_tool != tool:
                raise ValueError(
                    "Top-level 'ncrna_tool' differs between inputs:\n"
                    f"  - {tool_path or p}: {combined_tool}\n"
                    f"  - {p}: {tool}\n"
                    "This script requires ncrna_tool to be identical across input files."
                )

        features = cast(list[JsonValue], document["ncrna_features"])
        for raw_feat in features:
            n_features_in += 1
            feat = _coerce_ncrna_feature(raw_feat, source_path=p)
            lifted = _lift_feature_coords(
                feature=feat,
                chunk_re=chunk_re,
                agp_by_component=agp_by_component,
                allow_revcomp=allow_revcomp,
                source_path=p,
            )
            combined_features.append(lifted)

    if combined_analysis is None or combined_source is None or combined_tool is None:
        raise ValueError("No JSON files were read from the manifest (empty manifest?)")

    combined_json: dict[str, JsonValue] = {
        "analysis": cast(JsonValue, combined_analysis),
        "ncrna_tool": cast(JsonValue, combined_tool),
        "source": cast(JsonValue, combined_source),
        "ncrna_features": cast(JsonValue, combined_features),
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w", encoding="utf-8") as out_fh:
        json.dump(combined_json, out_fh, ensure_ascii=False, indent=2)
        out_fh.write("\n")

    schema_validator(out_json, json_schema=_FEATURE_SCHEMA_NAME)

    logging.info(f"Read {n_features_in} ncrna_feature(s); wrote {len(combined_features)} ncrna_feature(s)")


def combine_feature_json(
    json_manifest: Path,
    out_json: Path,
    chunk_re: re.Pattern[str],
    agp_file: Path | None = None,
    allow_revcomp: bool = False,
) -> None:
    """
    Combines JSON documents from split/chunked inputs and performs coordinate liftover.

    Supports:
    - load_repeat schema (repeat features + optional repeat consensus)
    - load_ncrna schema (ncRNA features for one tool type per file set)

    Liftover modes:
    1) AGP-driven (recommended):
        - Treat feature.seq_region as AGP component_id (part_id)
        - Lift seq_region_start/end into AGP object coordinates (record_start/end)
        - Replace seq_region with AGP object id (record)
        - If AGP orientation is '-', flip seq_region_strand (requires --allow-revcomp)

    2) Header-driven (no AGP):
        - Parse seq_region with --chunk-id-regex
        - Replace seq_region with base id and shift coordinates by (chunk_start - 1)

    Merging:
    - Features are concatenated after lifting.
    - Top-level analysis and source must be identical across input files.
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
        default=r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$",
        help=(
            "Regex used to identify chunked sequence IDs and extract coordinates when no AGP is provided. "
            "Must define named groups 'base' and 'start'."
        ),
    )

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
