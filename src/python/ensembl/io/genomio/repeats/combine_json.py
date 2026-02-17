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
"""Combine JSON payloads from split/chunked inputs and liftover coordinates, optionally using an AGP."""
import argparse
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
import gzip
import json
import logging
from pathlib import Path
import re
from typing import Any

from ensembl.io.genomio.schemas.json.validate import schema_validator
from ensembl.io.genomio.utils.agp_utils import AgpEntry, build_component_index, lift_range, parse_agp
from ensembl.io.genomio.utils.chunk_utils import get_paths_from_manifest, validate_regex

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


@dataclass(frozen=True)
class RecordLocation:
    path: Path
    description: str


def _pick_part_for_range(
    parts: list[AgpEntry], start: int, end: int, *, component_id: str, path: Path
) -> AgpEntry:
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


def _load_json_payload(path: Path) -> dict[str, Any]:
    """Load a single JSON object from (optionally gzipped) file."""
    with open_gz_file(path) as fh:
        raw = fh.read()

    if isinstance(raw, bytes):
        text = raw.decode("utf-8")
    else:
        text = raw

    obj = json.loads(text)
    if not isinstance(obj, dict):
        raise ValueError(f"Top-level JSON in {path} must be an object (got {type(obj).__name__})")
    return obj


def _schema_validate_file(path: Path, schema_name: str) -> None:
    """Validate a JSON file against a repository schema using schema_validator.

    schema_validator reads via Path(...).open(), so for .gz inputs we validate a temporary
    uncompressed copy in the same directory as output.
    """
    if path.suffix != ".gz":
        schema_validator(json_file=path, json_schema=schema_name)
        return

    # Decompress to a temp sibling file to keep schema_validator unchanged.
    tmp = path.with_suffix("")  # strip .gz
    # Avoid clobbering an existing real file.
    if tmp.exists():
        tmp = path.with_name(f"{tmp.name}.tmp_unzipped_for_schema_validation.json")

    logging.debug("Decompressing %s -> %s for schema validation", path, tmp)
    with gzip.open(path, "rt", encoding="utf-8") as in_fh, tmp.open("w", encoding="utf-8") as out_fh:
        out_fh.write(in_fh.read())

    try:
        schema_validator(json_file=tmp, json_schema=schema_name)
    finally:
        try:
            tmp.unlink(missing_ok=True)
        except Exception:
            logging.warning("Failed to remove temp validation file %s", tmp, exc_info=True)


def _lift_repeat_feature(
    feature: dict[str, Any],
    *,
    chunk_re: re.Pattern[str],
    agp_by_component: dict[str, list[AgpEntry]] | None,
    allow_revcomp: bool,
    source_path: Path,
) -> dict[str, Any]:
    """Return a modified copy of repeat_feature with lifted seq_region and coordinates."""
    seq_region = feature.get("seq_region")
    if not isinstance(seq_region, dict) or "name" not in seq_region:
        raise ValueError(f"repeat_feature.seq_region.name missing or invalid ({source_path})")

    name = seq_region["name"]
    if not isinstance(name, str) or not name:
        raise ValueError(f"repeat_feature.seq_region.name must be a non-empty string ({source_path})")

    sr_start = feature.get("seq_region_start")
    sr_end = feature.get("seq_region_end")
    strand = feature.get("seq_region_strand")

    if not isinstance(sr_start, int) or not isinstance(sr_end, int):
        raise ValueError(f"seq_region_start/end must be integers ({source_path})")
    if sr_start > sr_end:
        raise ValueError(f"seq_region_start > seq_region_end for {name}: {sr_start}>{sr_end} ({source_path})")
    if strand not in (-1, 1):
        raise ValueError(f"seq_region_strand must be -1 or 1 for {name} ({source_path})")

    out = dict(feature)
    out_seq_region = dict(seq_region)

    if agp_by_component is not None:
        parts = agp_by_component.get(name)
        if not parts:
            raise KeyError(f"seq_region '{name}' not found as an AGP component_id ({source_path})")

        part = _pick_part_for_range(parts, sr_start, sr_end, component_id=name, path=source_path)
        obj_id, new_start, new_end = lift_range(part, sr_start, sr_end, allow_revcomp)

        out_seq_region["name"] = obj_id
        out["seq_region"] = out_seq_region
        out["seq_region_start"] = new_start
        out["seq_region_end"] = new_end

        if part.orientation == "-":
            out["seq_region_strand"] = -strand

        return out

    # Header-driven (no AGP)
    m = chunk_re.match(name)
    if not m:
        return feature  # unchunked; nothing to do

    base = m.group("base")
    offset = int(m.group("start")) - 1

    out_seq_region["name"] = base
    out["seq_region"] = out_seq_region
    out["seq_region_start"] = sr_start + offset
    out["seq_region_end"] = sr_end + offset
    return out


def _merge_repeat_consensus(
    combined: dict[str, dict[str, Any]],
    incoming: list[dict[str, Any]],
    *,
    source_path: Path,
) -> None:
    """Merge repeat_consensus entries keyed by consensus_key; duplicates must be identical."""
    for rc in incoming:
        key = rc.get("consensus_key")
        if not isinstance(key, str) or not key:
            raise ValueError(f"repeat_consensus entry missing consensus_key ({source_path})")
        if key not in combined:
            combined[key] = rc
            continue
        if combined[key] != rc:
            raise ValueError(
                f"Conflicting repeat_consensus for consensus_key={key} ({source_path}). "
                "Entries differ; cannot safely merge."
            )


def _assert_same_top_level(
    meta_name: str, a: dict[str, Any], b: dict[str, Any], *, path_a: Path, path_b: Path
) -> None:
    if a != b:
        raise ValueError(
            f"Top-level '{meta_name}' differs between inputs:\n"
            f"  - {path_a}\n"
            f"  - {path_b}\n"
            "This script currently requires them to be identical."
        )


def combine_repeat_json(
    json_manifest: Path,
    out_json: Path,
    chunk_re: re.Pattern[str],
    json_schema: str,
    agp_file: Path | None = None,
    allow_revcomp: bool = False,
    validate_inputs: bool = True,
    validate_output: bool = True,
) -> None:
    """Combine split/chunked repeat JSON payloads back into original seq_region names and coordinates, optionally using an AGP.

    Coordinate lifting:
    1) AGP-driven (recommended):
        - Treat repeat_features[*].seq_region.name as AGP component_id (part_id)
        - Lift seq_region_start/end into AGP object coordinates (record_start/end)
        - Replace seq_region.name with AGP object id (record)
        - If AGP orientation is '-', flip seq_region_strand (requires --allow-revcomp)

    2) Header-driven (no AGP):
        - Parse seq_region.name with --chunk-id-regex
        - Replace seq_region.name with base id and shift coordinates by (chunk_start - 1)

    Merging:
    - repeat_consensus are merged by consensus_key; duplicates must be identical.
    - repeat_features are concatenated after lifting.
    - Top-level analysis and source must be identical across input files.
    """
    json_paths = get_paths_from_manifest(json_manifest)

    agp_by_component: dict[str, list[AgpEntry]] | None = None
    if agp_file is not None:
        agp_records = parse_agp(agp_file, allow_revcomp)
        agp_by_component = build_component_index(agp_records)

    combined_analysis: dict[str, Any] | None = None
    combined_source: dict[str, Any] | None = None
    analysis_path: Path | None = None
    source_path: Path | None = None

    consensus_by_key: dict[str, dict[str, Any]] = {}
    combined_features: list[dict[str, Any]] = []

    n_files = 0
    n_features_in = 0

    for p in json_paths:
        n_files += 1

        if validate_inputs:
            _schema_validate_file(p, json_schema)

        payload = _load_json_payload(p)

        analysis = payload["analysis"]
        source = payload["source"]

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

        _merge_repeat_consensus(consensus_by_key, payload["repeat_consensus"], source_path=p)

        for feat in payload["repeat_features"]:
            n_features_in += 1
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
        raise ValueError("No payloads were read from the manifest (empty manifest?)")

    # Ensure every feature's consensus_key exists in merged consensus
    missing: set[str] = set()
    for feat in combined_features:
        ck = feat["repeat_consensus_ref"]["consensus_key"]
        if ck not in consensus_by_key:
            missing.add(ck)
    if missing:
        missing_list = ", ".join(sorted(missing)[:20])
        extra = "" if len(missing) <= 20 else f" (+{len(missing) - 20} more)"
        raise ValueError(
            "repeat_features reference consensus_key(s) not present in repeat_consensus: "
            f"{missing_list}{extra}"
        )

    combined_payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "analysis": combined_analysis,
        "source": combined_source,
        "repeat_consensus": list(consensus_by_key.values()),
        "repeat_features": combined_features,
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w", encoding="utf-8") as out_fh:
        json.dump(combined_payload, out_fh, ensure_ascii=False, indent=2)
        out_fh.write("\n")

    if validate_output:
        schema_validator(json_file=out_json, json_schema=json_schema)

    logging.info("Read %d file(s)", n_files)
    logging.info(
        "Read %d repeat_feature(s); wrote %d repeat_feature(s)", n_features_in, len(combined_features)
    )
    logging.info("Wrote combined payload to %s", out_json)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--json-manifest",
        metavar="TXT",
        required=True,
        help="Manifest file containing paths of schema-valid JSON payloads to be combined (files may be gzipped).",
    )
    parser.add_argument_dst_path(
        "--out-json",
        metavar="JSON",
        required=True,
        help="Path for combined JSON payload output.",
    )
    parser.add_argument(
        "--json-schema",
        required=True,
        help=(
            "Schema identifier (as in ensembl.io.genomio.data.schemas) OR path to a JSON schema file, "
            "passed to schema_validator."
        ),
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
        "--no-validate-inputs",
        action="store_true",
        help="Skip schema validation for inputs (not recommended).",
    )
    parser.add_argument(
        "--no-validate-output",
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
            json_schema=args.json_schema,
            validate_inputs=not args.no_validate_inputs,
            validate_output=not args.no_validate_output,
        )
    except Exception:
        logging.exception(f"Error combining repeat JSON from files in {args.json_manifest}")
        raise


if __name__ == "__main__":
    main()
