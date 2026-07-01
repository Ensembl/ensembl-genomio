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
"""Parse Red repeat detector output into GenomIO repeat feature records."""

import argparse
from dataclasses import dataclass
from pathlib import Path

from ensembl.io.genomio.features.convert_to_genomio_json.base import (
    Consensus,
    format_parse_errors,
    has_valid_parsed_coordinates,
    parse_token,
)
from ensembl.io.genomio.features.convert_to_genomio_json.converters import (
    ConverterOptions,
    FeatureConverter,
    ParseFeaturesResult,
    _add_common_arguments,
)
from ensembl.utils.archive import open_gz_file

RED_RPT_COLUMNS = 3

RED_RPT_CONSENSUS = Consensus(
    name="repeatdetector",
    repeat_class="repeatdetector",
    repeat_type="repeatdetector",
    seq="N",
)

RED_RPT_CONSENSUS_KEY = RED_RPT_CONSENSUS.sha256_key()

__all__ = [
    "RedConverter",
    "RedParsedRow",
    "parse_red_data_row",
    "parse_red_output",
]


class RedConverter(FeatureConverter):
    """Converter for Red output."""

    analysis_logic_name = "repeatdetector"
    command = "red"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the Red subcommand parser."""
        red_parser = subparsers.add_parser(
            cls.command,
            help="Convert Red repeat detector output to GenomIO JSON.",
        )
        _add_common_arguments(red_parser)
        red_parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label="Repeats: Red",
            analysis_description=(
                'Repeats detected using <a href="https://bmcbioinformatics.biomedcentral.com/articles'
                '/10.1186/s12859-015-0654-5">Red (REPeatDetector)</a>'
            ),
            program="Red",
            repeatmasker_consensus_lib_path=None,
        )

    @classmethod
    def parse_features(
        cls,
        input_path: Path,
        _options: ConverterOptions | None = None,
    ) -> ParseFeaturesResult:
        """Parse Red output."""
        return parse_red_output(input_path)


@dataclass(frozen=True)
class RedParsedRow:
    """Parsed TRF row and its repeat consensus record."""

    feature: dict[str, object]


def parse_red_data_row(input_path: Path, line: str) -> RedParsedRow:
    """Parse a single TRF data row.

    Args:
        input_path: Input Red .rpt output path used for error messages.
        line: Raw Red data row without surrounding whitespace.

    Returns:
        Parsed row containing feature data and consensus.

    Raises:
        ValueError: If the row is malformed or contains invalid coordinates.

    """
    columns = line.split()
    if len(columns) != RED_RPT_COLUMNS:
        raise ValueError(
            f"Expected {RED_RPT_COLUMNS} columns in {input_path}, got {len(columns)}: line={line!r}"
        )

    seq_region = columns[0]
    seq_region_start = parse_token(int, columns[1], "start", line, input_path)
    seq_region_end = parse_token(int, columns[2], "end", line, input_path)
    repeat_length = (seq_region_end - seq_region_start) + 1

    has_valid_parsed_coordinates(
        input_path,
        seq_region_start=seq_region_start,
        seq_region_end=seq_region_end,
        repeat_start=1,
        repeat_end=repeat_length,
        line=line,
    )

    return RedParsedRow(
        feature={
            "seq_region": seq_region,
            "seq_region_start": seq_region_start,
            "seq_region_end": seq_region_end,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": repeat_length,
            "repeat_consensus": RED_RPT_CONSENSUS_KEY,
        }
    )


def parse_red_output(input_path: Path) -> ParseFeaturesResult:
    """Parse a Red .rpt file into repeat feature dictionaries and consensus records.

    Returns:
        A tuple containing:
            - repeat feature dictionaries
            - repeat consensus dictionary keyed by consensus SHA256 digest

    Raises:
        ValueError: If the TRF output contains malformed rows or invalid coordinate values.

    """
    features: list[dict] = []
    consensuses_by_key: dict[str, Consensus] = {RED_RPT_CONSENSUS_KEY: RED_RPT_CONSENSUS}
    errors: list[str] = []

    with open_gz_file(input_path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue

            try:
                parsed_row = parse_red_data_row(input_path, line)
            except ValueError as exc:
                errors.append(str(exc))
                continue

            features.append(parsed_row.feature)

    if errors:
        raise ValueError(format_parse_errors("Red output", input_path, errors))

    return features, consensuses_by_key
