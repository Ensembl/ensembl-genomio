#!/usr/bin/env python3

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Parse tRNAscan output into GenomIO JSON converter records."""

__all__ = [
    "TrnaScanConverter",
]

import argparse
from datetime import datetime, timezone
import json
import re
import sys
from pathlib import Path
from typing import Any, Optional

BASE_COMPLEMENTS = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def anticodon_to_codon(anticodon: str) -> str:
    """Return the DNA codon paired with a tRNA anticodon."""
    return anticodon.translate(BASE_COMPLEMENTS)[::-1].upper()


def is_pseudogene(trna_type: str, note: Optional[str]) -> bool:
    """Detect pseudogene annotations from tRNAscan-SE fields."""
    values = [trna_type]
    if note:
        values.append(note)
    return any("pseudo" in value.lower() for value in values)


def parse_trnascan_line(line: str, line_number: int) -> Optional[dict[str, Any]]:
    """Parse one tRNAscan-SE output row."""
    stripped = line.strip()
    if not stripped or stripped.startswith("-"):
        return None

    fields = [field.strip() for field in (stripped.split("\t") if "\t" in stripped else re.split(r"\s+", stripped))]
    if len(fields) < 9:
        return None
    if not fields[1].isdigit() or not fields[2].isdigit() or not fields[3].isdigit():
        return None

    try:
        begin = int(fields[2])
        end = int(fields[3])
        isotype = fields[4]
        anticodon = fields[5]
        note = " ".join(fields[9:]) if len(fields) > 9 else None

        return {
            "seq_region": fields[0],
            "seq_region_start": min(begin, end),
            "seq_region_end": max(begin, end),
            "seq_region_strand": "+" if begin <= end else "-",
            "biotype": "tRNA",
            "display_label": f"tRNA-{isotype}",
            "score": float(fields[8]),
            "isotype": isotype,
            "anticodon": anticodon,
            "codon": anticodon_to_codon(anticodon),
            "is_pseudogene": is_pseudogene(isotype, note),
        }
    except ValueError as exc:
        raise ValueError(f"Could not parse data row on line {line_number}: {line.rstrip()}") from exc


def parse_trnascan_features(input_tsv: Path) -> list[dict[str, Any]]:
    """Parse tRNAscan-SE tabular output into feature records."""
    records: list[dict[str, Any]] = []

    with input_tsv.open("r", encoding="utf-8") as in_handle:
        for line_number, line in enumerate(in_handle, start=1):
            record = parse_trnascan_line(line, line_number)
            if record is not None:
                records.append(record)

    return records


def _add_trnascan_arguments(subparser: argparse.ArgumentParser) -> None:
    """Add tRNAscan-SE arguments, with a plain-argparse fallback for local execution."""
    try:
        subparser.add_argument_src_path("--input", required=True, help="Input tRNAscan-SE TSV file")
        subparser.add_argument_dst_path("--output", required=True, help="JSON output path")
    except AttributeError:
        subparser.add_argument("--input", required=True, type=Path, help="Input tRNAscan-SE TSV file")
        subparser.add_argument("--output", required=True, type=Path, help="JSON output path")


class TrnaScanConverter:
    """Converter for tRNAscan-SE output."""

    analysis_logic_name = "trnascan"
    analysis_display_label = "tRNAs"
    analysis_description = "tRNA genes predicted by tRNAscan-SE."
    command = "trnascan"
    program = "tRNAscan-SE"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the tRNAscan-SE subcommand parser."""
        trnascan_parser = subparsers.add_parser(
            cls.command,
            help="Convert tRNAscan-SE output to GenomIO JSON.",
        )
        _add_trnascan_arguments(trnascan_parser)
        trnascan_parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label=cls.analysis_display_label,
            analysis_description=cls.analysis_description,
            program=cls.program,
        )

    @classmethod
    def parse_features(
        cls,
        input_path: Path,
        _options: Optional[object] = None,
    ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
        """Parse tRNAscan-SE output."""
        return parse_trnascan_features(input_path), {}


def build_genomio_json_file(input_tsv: Path) -> dict[str, object]:
    """Build a GenomIO JSON file from tRNAscan-SE results."""
    return {
        "analysis": {
            "run_date": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "logic_name": "trnascan",
            "display_label": "tRNAs",
            "description": "tRNA genes predicted by tRNAscan-SE.",
            "program": "tRNAscan-SE",
        },
        "source": {
            "source_provider": "tRNAscan-SE",
            "is_primary": True,
        },
        "ncrna_tool": "trnascan",
        "ncrna_features": parse_trnascan_features(input_tsv),
    }


def write_trnascan_json(input_tsv: Path, output_json: Path) -> None:
    """Read tRNAscan-SE tabular output and write the GenomIO JSON file."""
    document = build_genomio_json_file(input_tsv)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def add_trnascan_arguments(subparser: argparse.ArgumentParser) -> None:
    """Add tRNAscan-SE specific CLI arguments."""
    subparser.add_argument("--input", required=True, type=Path, help="Input tRNAscan-SE TSV file")
    subparser.add_argument("--output", required=True, type=Path, help="JSON output path")


def build_parser() -> argparse.ArgumentParser:
    """Build the tRNAscan-SE converter CLI parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="tool", required=True)
    trnascan_parser = subparsers.add_parser("trnascan", help="Convert tRNAscan-SE output to GenomIO JSON.")
    add_trnascan_arguments(trnascan_parser)
    trnascan_parser.set_defaults(handler=run_trnascan)
    return parser


def parse_args(arg_list: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    if arg_list is None:
        arg_list = sys.argv[1:]
    if arg_list is not None and arg_list and arg_list[0] != "trnascan":
        arg_list = ["trnascan", *arg_list]
    return build_parser().parse_args(arg_list)


def run_trnascan(args: argparse.Namespace) -> int:
    """Run the tRNAscan-SE conversion command."""
    try:
        write_trnascan_json(input_tsv=args.input, output_json=args.output)
    except Exception as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


def main(arg_list: Optional[list[str]] = None) -> int:
    """Command-line entry point."""
    args = parse_args(arg_list)
    handler = getattr(args, "handler", None)
    if handler is None:
        raise ValueError("No command selected")
    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
