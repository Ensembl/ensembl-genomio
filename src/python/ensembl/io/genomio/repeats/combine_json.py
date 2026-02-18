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

"""Combines JSON documents from split/chunked inputs and performs coordinate liftover, optionally using an AGP."""

import argparse
from collections.abc import Iterable
import json
import logging
from pathlib import Path
import re
from typing import NotRequired, TypeAlias, TypedDict, cast

from ensembl.io.genomio.schemas.json.validate import schema_validator
from ensembl.io.genomio.utils.agp_utils import AgpEntry, build_component_index, lift_range, parse_agp
from ensembl.io.genomio.utils.chunk_utils import get_paths_from_manifest, validate_regex

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


_REPEAT_SCHEMA_NAME = "repeat"

JsonValue: TypeAlias = None | bool | int | float | str | list["JsonValue"] | dict[str, "JsonValue"]


class RepeatConsensus(TypedDict):
    repeat_name: str
    repeat_class: str
    repeat_type: str
    repeat_consensus: NotRequired[str]


class RepeatFeature(TypedDict):
    seq_region: str
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: int  # constrained to -1/+1 by validation
    repeat_start: int
    repeat_end: int
    repeat_consensus: NotRequired[RepeatConsensus | None]
    score: NotRequired[float | None]
    attributes: NotRequired[dict[str, JsonValue]]


def _get_agp_entry_for_range(
    parts: list[AgpEntry], start: int, end: int, *, component_id: str, path: Path
) -> AgpEntry:
    """Retrieves the unique AGP part whose component span fully contains a given range."""
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
    if not isinstance(parsed, dict):
        raise ValueError(f"Top-level JSON in {path} must be an object (got {type(parsed).__name__})")

    return cast(dict[str, JsonValue], parsed)


def _schema_validate_file(path: Path) -> None:
    """Validates a JSON file against a repository schema using schema_validator."""
    if path.suffix != ".gz":
        schema_validator(json_file=path, json_schema=_REPEAT_SCHEMA_NAME)
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
        schema_validator(json_file=tmp, json_schema=_REPEAT_SCHEMA_NAME)
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


def _coerce_repeat_consensus(obj: JsonValue, source_path: Path) -> RepeatConsensus:
    """Validates and casts an arbitrary object to RepeatConsensus."""
    if not isinstance(obj, dict):
        raise ValueError(f"repeat_consensus entries must be objects ({source_path})")

    rn = obj.get("repeat_name")
    rc_class = obj.get("repeat_class")
    rt = obj.get("repeat_type")

    if not isinstance(rn, str) or not rn:
        raise ValueError(f"repeat_consensus.repeat_name must be a non-empty string ({source_path})")
    if not isinstance(rc_class, str) or not rc_class:
        raise ValueError(f"repeat_consensus.repeat_class must be a non-empty string ({source_path})")
    if not isinstance(rt, str) or not rt:
        raise ValueError(f"repeat_consensus.repeat_type must be a non-empty string ({source_path})")

    if "repeat_consensus" in obj:
        rc_seq = obj.get("repeat_consensus")
        if rc_seq is not None and not isinstance(rc_seq, str):
            raise ValueError(f"repeat_consensus.repeat_consensus must be a string if present ({source_path})")

    return cast(RepeatConsensus, obj)


def _coerce_repeat_feature(obj: JsonValue, source_path: Path) -> RepeatFeature:
    """Validates and casts an arbitrary object to RepeatFeature (minimal fields)."""
    if not isinstance(obj, dict):
        raise ValueError(f"repeat_features entries must be objects ({source_path})")

    seq_region = obj.get("seq_region")
    if not isinstance(seq_region, str) or not seq_region:
        raise ValueError(f"repeat_feature.seq_region must be a non-empty string ({source_path})")

    for k in ("seq_region_start", "seq_region_end", "repeat_start", "repeat_end"):
        v = obj.get(k)
        if not isinstance(v, int):
            raise ValueError(f"repeat_feature.{k} must be an integer ({source_path})")

    strand = obj.get("seq_region_strand")
    if strand not in (-1, 1):
        raise ValueError(f"repeat_feature.seq_region_strand must be -1 or 1 ({source_path})")

    if "score" in obj:
        score = obj.get("score")
        if score is not None and not isinstance(score, (int, float)):
            raise ValueError(f"repeat_feature.score must be numeric or null ({source_path})")

    if "attributes" in obj:
        attrs = obj.get("attributes")
        if attrs is not None and not isinstance(attrs, dict):
            raise ValueError(f"repeat_feature.attributes must be an object if present ({source_path})")

    if "repeat_consensus" in obj:
        rc_ref = obj.get("repeat_consensus")
        if rc_ref is not None:
            _coerce_repeat_consensus(rc_ref, source_path=source_path)

    return cast(RepeatFeature, obj)


def _consensus_key(rc: RepeatConsensus, source_path: Path) -> tuple[str, str, str]:
    """Extracts and validates the natural key for a repeat_consensus entry."""
    rn = rc["repeat_name"]
    rc_class = rc["repeat_class"]
    rt = rc["repeat_type"]

    if not rn or not rc_class or not rt:
        raise ValueError(
            "repeat_consensus entries must include non-empty repeat_name, repeat_class, repeat_type "
            f"({source_path})"
        )
    return rn, rc_class, rt


def _merge_repeat_consensus(
    combined: dict[tuple[str, str, str], RepeatConsensus],
    incoming: Iterable[JsonValue],
    source_path: Path,
) -> None:
    """Merges repeat_consensus entries keyed by (repeat_name, repeat_class, repeat_type); duplicates must be identical."""
    for raw in incoming:
        rc = _coerce_repeat_consensus(raw, source_path=source_path)
        key = _consensus_key(rc, source_path=source_path)
        if key not in combined:
            combined[key] = rc
            continue
        if combined[key] != rc:
            raise ValueError(
                f"Conflicting repeat_consensus for key={key} ({source_path}). "
                "Entries differ; cannot safely merge."
            )


def _feature_consensus_key(feature: RepeatFeature, source_path: Path) -> tuple[str, str, str] | None:
    """Constructs the natural consensus key from a feature."""
    ref = feature.get("repeat_consensus")
    if ref is None:
        return None

    rn = ref["repeat_name"]
    rc_class = ref["repeat_class"]
    rt = ref["repeat_type"]

    if not rn or not rc_class or not rt:
        raise ValueError(
            "repeat_feature.repeat_consensus must contain non-empty repeat_name, repeat_class, repeat_type "
            f"({source_path})"
        )

    return rn, rc_class, rt


def _lift_repeat_feature(
    feature: RepeatFeature,
    chunk_re: re.Pattern[str],
    agp_by_component: dict[str, list[AgpEntry]] | None,
    allow_revcomp: bool,
    source_path: Path,
) -> RepeatFeature:
    """Constructs a modified copy of a repeat_feature with lifted seq_region and coordinates."""
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

        return cast(RepeatFeature, out)

    m = chunk_re.match(name)
    if not m:
        return feature  # unchunked; nothing to do

    base = m.group("base")
    offset = int(m.group("start")) - 1

    out["seq_region"] = base
    out["seq_region_start"] = sr_start + offset
    out["seq_region_end"] = sr_end + offset
    return cast(RepeatFeature, out)


def combine_repeat_json(
    json_manifest: Path,
    out_json: Path,
    chunk_re: re.Pattern[str],
    agp_file: Path | None = None,
    allow_revcomp: bool = False,
    validate_inputs: bool = True,
    validate_output: bool = True,
) -> None:
    """
    Combines JSON documents from split/chunked inputs and performs coordinate liftover.

    JSON inputs must conform to the schema defined in ensembl.io.genomio.data.schemas.repeat.

    Repeat feature lifting:
    1) AGP-driven (recommended):
        - Treat repeat_features[*].seq_region as AGP component_id (part_id)
        - Lift seq_region_start/end into AGP object coordinates (record_start/end)
        - Replace seq_region with AGP object id (record)
        - If AGP orientation is '-', flip seq_region_strand (requires --allow-revcomp)

    2) Header-driven (no AGP):
        - Parse seq_region with --chunk-id-regex
        - Replace seq_region with base id and shift coordinates by (chunk_start - 1)

    Merging:
    - repeat_consensus are merged by natural key (repeat_name, repeat_class, repeat_type); duplicates must be identical.
    - repeat_features are concatenated after lifting.
    - Top-level analysis and source must be identical across input files.
    - If a feature has a non-null repeat_consensus reference, it must exist in the merged repeat_consensus list.
    """
    json_paths = get_paths_from_manifest(json_manifest)

    agp_by_component: dict[str, list[AgpEntry]] | None = None
    if agp_file is not None:
        agp_records = parse_agp(agp_file, allow_revcomp)
        agp_by_component = build_component_index(agp_records)

    combined_analysis: dict[str, JsonValue] | None = None
    combined_source: dict[str, JsonValue] | None = None
    analysis_path: Path | None = None
    source_path: Path | None = None

    consensus_by_key: dict[tuple[str, str, str], RepeatConsensus] = {}
    combined_features: list[RepeatFeature] = []

    n_files = 0
    n_features_in = 0

    for p in json_paths:
        n_files += 1

        if validate_inputs:
            _schema_validate_file(p)

        document = _load_json_document(p)

        analysis_val = document.get("analysis")
        source_val = document.get("source")

        if not isinstance(analysis_val, dict):
            raise ValueError(f"Top-level 'analysis' must be an object ({p})")
        if not isinstance(source_val, dict):
            raise ValueError(f"Top-level 'source' must be an object ({p})")

        analysis = cast(dict[str, JsonValue], analysis_val)
        source = cast(dict[str, JsonValue], source_val)

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

        incoming_consensus_val = document.get("repeat_consensus", [])
        if incoming_consensus_val is None:
            incoming_consensus_val = []
        if not isinstance(incoming_consensus_val, list):
            raise ValueError(f"Top-level repeat_consensus must be a list if present ({p})")

        _merge_repeat_consensus(consensus_by_key, incoming_consensus_val, source_path=p)

        feats_val = document.get("repeat_features")
        if not isinstance(feats_val, list):
            raise ValueError(f"Top-level repeat_features must be a list ({p})")

        for raw_feat in feats_val:
            n_features_in += 1
            feat = _coerce_repeat_feature(raw_feat, source_path=p)
            combined_features.append(
                _lift_repeat_feature(
                    feat,
                    chunk_re=chunk_re,
                    agp_by_component=agp_by_component,
                    allow_revcomp=allow_revcomp,
                    source_path=p,
                )
            )

    if combined_analysis is None or combined_source is None:
        raise ValueError("No JSON files were read from the manifest (empty manifest?)")

    missing: set[tuple[str, str, str]] = set()
    for feat in combined_features:
        key = _feature_consensus_key(feat, source_path=out_json)
        if key is None:
            continue
        if key not in consensus_by_key:
            missing.add(key)

    if missing:
        missing_list = ", ".join([f"{k[0]}|{k[1]}|{k[2]}" for k in sorted(missing)[:20]])
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

    if validate_output:
        schema_validator(out_json, json_schema=_REPEAT_SCHEMA_NAME)

    logging.info(f"Read {n_files} file(s)")
    logging.info(f"Read {n_features_in} repeat_feature(s); wrote {len(combined_features)} repeat_feature(s)")
    logging.info(f"Wrote combined file to {out_json}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--json-manifest",
        metavar="TXT",
        required=True,
        help="Manifest file containing paths of schema-valid JSON files to be combined (files may be gzipped).",
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
        default=None,
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
    parser.add_argument(
        "--no-input-validation",
        action="store_true",
        help="Skip schema validation for inputs (not recommended).",
    )
    parser.add_argument(
        "--no-output-validation",
        action="store_true",
        help="Skip schema validation for output (not recommended).",
    )

    args = parser.parse_args(argv)
    init_logging_with_args(args)
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    chunk_re = validate_regex(args.chunk_id_regex)
    try:
        combine_repeat_json(
            json_manifest=args.json_manifest,
            out_json=args.out_json,
            agp_file=args.agp_file,
            allow_revcomp=args.allow_revcomp,
            chunk_re=chunk_re,
            validate_inputs=not args.no_input_validation,
            validate_output=not args.no_output_validation,
        )
    except Exception:
        logging.exception(f"Error combining repeat JSON from files in {args.json_manifest}")
        raise


if __name__ == "__main__":
    main()
