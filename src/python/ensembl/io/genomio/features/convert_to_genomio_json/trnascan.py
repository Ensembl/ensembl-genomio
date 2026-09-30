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
"""Parse tRNAscan-SE output into GenomIO ncRNA feature records."""

__all__ = [
    "TrnaScanConverter",
    "TrnaScanParsedRow",
]

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

from ensembl.io.genomio.features.convert_to_genomio_json.base import (
    ConverterOptions,
    FeatureConverter,
    ParseFeaturesResult,
    format_parse_errors,
    parse_token,
    register_converter,
    register_top_level_converter,
)

TRNASCAN_MIN_COLUMNS = 9
TRNASCAN_HEADER_NAMES = {"Sequence", "Name"}
BASE_COMPLEMENTS = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


@register_top_level_converter
@register_converter
class TrnaScanConverter(FeatureConverter):
    """Converter for tRNAscan-SE output."""

    analysis_logic_name = "trnascan"
    command = "trnascan"
    feature_collection_name = "ncrna_features"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the tRNAscan-SE subcommand parser."""
        trnascan_parser = subparsers.add_parser(
            cls.command,
            help="Convert tRNAscan-SE output to GenomIO JSON.",
        )
        cls.add_common_arguments(trnascan_parser)
        trnascan_parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label="tRNAs",
            analysis_description="tRNA genes predicted by tRNAscan-SE.",
            program="tRNAscan-SE",
        )

    @classmethod
    def parse_features(
        cls,
        input_path: Path,
        _options: ConverterOptions | None = None,
    ) -> ParseFeaturesResult:
        """Parse tRNAscan-SE output."""
        return parse_output(input_path)

    @classmethod
    def additional_json_fields(cls) -> dict[str, object]:
        """Add tRNAscan-specific metadata to the JSON document."""
        return {"ncrna_tool": "trnascan"}

@dataclass(frozen=True)
class TrnaScanParsedRow:
    """Parsed tRNAscan-SE row."""

    feature: dict[str, object]


def anticodon_to_codon(anticodon: str) -> str:
    """Return the DNA codon paired with a tRNA anticodon."""
    return anticodon.translate(BASE_COMPLEMENTS)[::-1].upper()


def is_pseudogene(trna_type: str, note: str | None) -> bool:
    """Detect pseudogene annotations from tRNAscan-SE fields."""
    values = [trna_type]
    if note:
        values.append(note)
    return any("pseudo" in value.lower() for value in values)


def parse_row(input_path: Path, line: str) -> TrnaScanParsedRow:
    """Parse a single tRNAscan-SE data row.

    Args:
        input_path: Input path used in parsing error messages.
        line: Raw tRNAscan-SE data row without surrounding whitespace.

    Returns:
        Parsed row containing one ncRNA feature.

    Raises:
        ValueError: If the row is malformed or contains invalid numeric values.

    """
    columns = [column.strip() for column in (line.split("\t") if "\t" in line else re.split(r"\s+", line))]
    if len(columns) < TRNASCAN_MIN_COLUMNS:
        raise ValueError(
            f"Expected at least {TRNASCAN_MIN_COLUMNS} columns in {input_path}, "
            f"got {len(columns)}: line={line!r}"
        )

    parse_token(int, columns[1], "tRNA number", line, input_path)
    begin = parse_token(int, columns[2], "begin coordinate", line, input_path)
    end = parse_token(int, columns[3], "end coordinate", line, input_path)
    score = parse_token(float, columns[8], "score", line, input_path)

    if begin < 1 or end < 1:
        raise ValueError(
            f"Invalid coordinates in {input_path}: begin={begin}, end={end}, line={line!r}"
        )

    isotype = columns[4]
    anticodon = columns[5]
    note = " ".join(columns[TRNASCAN_MIN_COLUMNS:]) or None

    return TrnaScanParsedRow(
        feature={
            "seq_region": columns[0],
            "seq_region_start": min(begin, end),
            "seq_region_end": max(begin, end),
            "seq_region_strand": "+" if begin <= end else "-",
            "biotype": "tRNA",
            "display_label": f"tRNA-{isotype}",
            "score": score,
            "isotype": isotype,
            "anticodon": anticodon,
            "codon": anticodon_to_codon(anticodon),
            "is_pseudogene": is_pseudogene(isotype, note),
        }
    )


def parse_output(input_path: Path) -> ParseFeaturesResult:
    """Parse a tRNAscan-SE output file into ncRNA feature dictionaries.

    Returns:
        A tuple containing the ncRNA features and an empty consensus mapping.

    Raises:
        ValueError: If one or more data rows are malformed.

    """
    features: list[dict[str, object]] = []
    errors: list[str] = []

    with input_path.open("r", encoding="utf-8") as input_handle:
        for raw_line in input_handle:
            line = raw_line.strip()
            if not line or line.startswith("-"):
                continue

            first_column = line.split("\t", maxsplit=1)[0].split(maxsplit=1)[0]
            if first_column in TRNASCAN_HEADER_NAMES:
                continue

            try:
                parsed_row = parse_row(input_path, line)
            except ValueError as exc:
                errors.append(str(exc))
                continue

            features.append(parsed_row.feature)

    if errors:
        raise ValueError(format_parse_errors("tRNAscan-SE output", input_path, errors))

    return features, {}

