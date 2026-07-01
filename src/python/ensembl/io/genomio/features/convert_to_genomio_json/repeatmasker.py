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
"""Parse RepeatMasker output into GenomIO repeat feature records."""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

from Bio import SeqIO

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
from ensembl.utils.argparse import ArgumentParser

REPEATMASKER_MAPPINGS = [
    (r"^Low_Comp.*$", "Low complexity regions"),
    (r"^LINE.*$", "Type I Transposons/LINE"),
    (r"^SINE.*$", "Type I Transposons/SINE"),
    (r"^DNA.*$", "Type II Transposons"),
    (r"^LTR.*$", "LTRs"),
    (r"^Other.*$", "Other repeats"),
    (r"^Satelli.*$", "Satellite repeats"),
    (r"^Simple.*$", "Simple repeats"),
    (r"^Tandem.*$", "Tandem repeats"),
    (r"^TRF.*$", "Tandem repeats"),
    (r"^Waterman$", "Waterman"),
    (r"^Recon$", "Recon"),
    (r"^Tet_repeat$", "Tetraodon repeats"),
    (r"^MaskRegion$", "Mask region"),
    (r"^dust.*$", "Dust"),
    (r"^Unknown.*$", "Unknown"),
    (r".*RNA$", "RNA repeats"),
]

REPEATMASKER_COMPILED_MAPPINGS = [(re.compile(pattern), mapped) for pattern, mapped in REPEATMASKER_MAPPINGS]

__all__ = [
    "RepeatMaskerConverter",
    "RepeatMaskerCustomConverter",
    "RepeatMaskerParsedRow",
    "RepeatMaskerRepbaseConverter",
]


def _add_repeatmasker_common_args(subparser: ArgumentParser) -> None:
    """Add arguments shared by RepeatMasker subcommands."""
    _add_common_arguments(subparser)
    subparser.add_argument_src_path(
        "--consensus-lib",
        metavar="RM_LIB",
        default=argparse.SUPPRESS,
        help="FASTA file containing consensus sequences for the RepeatMasker library used.",
    )


class RepeatMaskerCustomConverter(FeatureConverter):
    """Converter for RepeatMasker output generated with a custom library."""

    analysis_logic_name = "repeatmask_customlib"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the custom-library RepeatMasker mode parser."""
        custom_parser = subparsers.add_parser(
            "custom",
            help="Convert RepeatMasker output generated using a custom library.",
        )
        _add_repeatmasker_common_args(custom_parser)
        custom_parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label="Repeats: Custom library",
            analysis_description=(
                'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>,'
                " using a custom library of <em>ab initio</em> repeat profiles for this species."
            ),
        )

    @classmethod
    def parse_features(
        cls,
        input_path: Path,
        options: ConverterOptions | None = None,
    ) -> ParseFeaturesResult:
        """Parse custom-library RepeatMasker output."""
        options = options or ConverterOptions()
        return parse_output(input_path, options.repeatmasker_consensus_lib_path)


class RepeatMaskerRepbaseConverter(RepeatMaskerCustomConverter):
    """Converter for RepeatMasker output generated with Repbase."""

    analysis_logic_name = "repeatmask_repbase"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the Repbase RepeatMasker mode parser."""
        repbase_parser = subparsers.add_parser(
            "repbase",
            help="Convert RepeatMasker output geneterated using Repbase.",
        )
        _add_repeatmasker_common_args(repbase_parser)
        repbase_parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label="Repeats: Repbase",
            analysis_description=(
                'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>,'
                ' using the <a rel="external" href="http://www.girinst.org/repbase/">Repbase</a> library of '
                "repeat profiles."
            ),
        )


class RepeatMaskerConverter:
    """Top-level RepeatMasker CLI converter."""

    command = "repeatmasker"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the RepeatMasker subcommand parser."""
        repeatmasker_parser = subparsers.add_parser(
            cls.command,
            help="Convert RepeatMasker output to GenomIO JSON.",
        )
        repeatmasker_parser.set_defaults(program="RepeatMasker")
        repeatmasker_subparsers = repeatmasker_parser.add_subparsers(dest="repeatmasker_mode", required=True)
        RepeatMaskerCustomConverter.add_parser(repeatmasker_subparsers)
        RepeatMaskerRepbaseConverter.add_parser(repeatmasker_subparsers)


@dataclass(frozen=True)
class RepeatMaskerParsedRow:
    """Parsed RepeatMasker row and its consensus lookup triplet."""

    feature: dict[str, object]
    consensus_triplet: tuple[str, str, str]


def map_repeatmasker_repeat_consensus_type(repeat_class: str) -> str:
    """Map a raw RepeatMasker repeat class to a GenomIO repeat category.

    Args:
        repeat_class: Raw repeat class string extracted from RepeatMasker output.

    Returns:
        Mapped repeat type category. Returns "Unknown" when no mapping matches.

    """
    for regex, mapped in REPEATMASKER_COMPILED_MAPPINGS:
        if regex.match(repeat_class):
            return mapped
    return "Unknown"


def parse_repeatmasker_consensus_library(
    consensus_lib_path: Path,
) -> tuple[dict[tuple[str, str, str], str], dict[str, Consensus]]:
    """Parse a RepeatMasker consensus library FASTA file into a dictionary of Consensus records.

    The parser expects FASTA headers in the format:
        >consensus_name#repeat_class/repeat_type

    Args:
        consensus_lib_path: Path to the RepeatMasker consensus library FASTA file.

    Returns:
        A tuple containing:
            - A dictionary mapping (consensus_name, repeat_class, repeat_type) tuples to SHA256 digests.
            - A dictionary of Consensus records keyed by SHA256 digest.

    """
    consensus_keys_by_triplet: dict[tuple[str, str, str], str] = {}
    consensuses_by_key: dict[str, Consensus] = {}
    with open_gz_file(consensus_lib_path) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            repeat_name = record.id
            repeat_class = "Unknown"
            repeat_type = "Unknown"

            if "#" in record.id:
                repeat_name, class_type_part = record.id.split("#", 1)
                if "/" in class_type_part:
                    repeat_class, repeat_type = class_type_part.split("/", 1)
                else:
                    repeat_class = class_type_part

            consensus_obj = Consensus(
                name=repeat_name,
                repeat_class=repeat_class,
                repeat_type=map_repeatmasker_repeat_consensus_type(repeat_class),
                seq=str(record.seq),
            )
            consensus_key = consensus_obj.sha256_key()
            consensus_keys_by_triplet[(repeat_name, repeat_class, repeat_type)] = consensus_key
            consensuses_by_key[consensus_key] = consensus_obj
    return consensus_keys_by_triplet, consensuses_by_key


def parse_repeat_class_field(
    input_path: Path,
    repeat_class_field: str,
    line: str,
) -> tuple[str, str]:
    """Parse a RepeatMasker repeat class/family field into class and raw repeat type.

    Args:
        input_path: Input RepeatMasker output path used for error messages.
        repeat_class_field: Raw class/family field from a RepeatMasker row.
        line: Original input line.

    Returns:
        Repeat class and raw repeat type.

    Raises:
        ValueError: If a slash-delimited class/family field is malformed.

    """
    if "/" not in repeat_class_field:
        return repeat_class_field, "Unknown"

    repeat_class, repeat_type = repeat_class_field.split("/", 1)
    if not repeat_class or not repeat_type:
        raise ValueError(
            f"Malformed repeat_class/family in {input_path}: value={repeat_class_field}, line={line!r}"
        )
    return repeat_class, repeat_type


def parse_strand_coordinates(
    input_path: Path,
    columns: list[str],
    line: str,
) -> tuple[str, int, int]:
    """Parse RepeatMasker strand and repeat coordinates.

    Args:
        input_path: Input RepeatMasker output path used for error messages.
        columns: Split RepeatMasker row columns.
        line: Original input line.

    Returns:
        Sequence-region strand, repeat start, and repeat end.

    Raises:
        ValueError: If the strand token or coordinate tokens are invalid.

    """
    strand_token = columns[8]
    if strand_token == "+":
        return (
            "+",
            parse_token(int, columns[11], "repeat_start", line, input_path),
            parse_token(int, columns[12], "repeat_end", line, input_path),
        )
    if strand_token == "C":
        return (
            "-",
            parse_token(int, columns[13], "repeat_start", line, input_path),
            parse_token(int, columns[12], "repeat_end", line, input_path),
        )
    raise ValueError(f"Unexpected strand token in {input_path}: token={strand_token!r}, line={line!r}")


def parse_row(input_path: Path, line: str) -> RepeatMaskerParsedRow | None:
    """Parse a single RepeatMasker data row.

    Args:
        input_path: Input RepeatMasker output path used for error messages.
        line: Raw RepeatMasker row without surrounding whitespace.

    Returns:
        Parsed row, or ``None`` if repeat coordinates are invalid and the row should be skipped.

    Raises:
        ValueError: If the row is malformed or contains invalid sequence-region coordinates.

    """
    columns = line.split()
    if columns[-1] == "*":
        columns.pop()

    if len(columns) < 14 or len(columns) > 15:  # noqa: PLR2004  -- ignore ruff "magic value" rule
        raise ValueError(f"Expected 14 or 15 columns in {input_path}, got {len(columns)}: line={line!r}")

    score = parse_token(float, columns[0], "score", line, input_path)
    perc_div = parse_token(float, columns[1], "perc_div", line, input_path)
    perc_del = parse_token(float, columns[2], "perc_del", line, input_path)
    perc_ins = parse_token(float, columns[3], "perc_ins", line, input_path)

    seq_region = columns[4]
    seq_region_start = parse_token(int, columns[5], "seq_region_start", line, input_path)
    seq_region_end = parse_token(int, columns[6], "seq_region_end", line, input_path)
    seq_region_strand, repeat_start, repeat_end = parse_strand_coordinates(
        input_path,
        columns,
        line,
    )
    repeat_name = columns[9]
    repeat_class, repeat_type = parse_repeat_class_field(input_path, columns[10], line)

    if not has_valid_parsed_coordinates(
        input_path,
        seq_region_start=seq_region_start,
        seq_region_end=seq_region_end,
        repeat_start=repeat_start,
        repeat_end=repeat_end,
        line=line,
    ):
        return None

    return RepeatMaskerParsedRow(
        feature={
            "seq_region": seq_region,
            "seq_region_start": seq_region_start,
            "seq_region_end": seq_region_end,
            "seq_region_strand": seq_region_strand,
            "repeat_start": repeat_start,
            "repeat_end": repeat_end,
            "score": score,
            "attributes": {
                "perc_div": perc_div,
                "perc_del": perc_del,
                "perc_ins": perc_ins,
                "repeatmasker_repeat_type": repeat_type,
            },
        },
        consensus_triplet=(repeat_name, repeat_class, repeat_type),
    )


def parse_output(input_path: Path, consensus_lib_path: Path | None) -> ParseFeaturesResult:
    """Parse a RepeatMasker .out file into repeat feature dictionaries and consensus records.

    If a RepeatMasker consensus library FASTA file is provided, consensus sequences will be
    extracted from the library and associated with features based on their repeat names.
    Otherwise, an empty consensus dictionary will be returned.

    Args:
        input_path: Path to the RepeatMasker .out file to parse.
        consensus_lib_path: Optional path to a RepeatMasker consensus library FASTA file.

    Returns:
        A tuple containing:
            - repeat feature dictionaries
            - repeat consensus dictionary keyed by consensus SHA256 digest

    Raises:
        ValueError: If the RepeatMasker output contains malformed rows or invalid coordinate values.

    """
    consensus_keys_by_triplet: dict[tuple[str, str, str], str] = {}
    consensuses_by_key: dict[str, Consensus] = {}
    if consensus_lib_path is not None:
        consensus_keys_by_triplet, consensuses_by_key = parse_repeatmasker_consensus_library(
            consensus_lib_path
        )

    features: list[dict] = []
    errors: list[str] = []
    with open_gz_file(input_path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue

            lower_line = line.lower()
            if (
                "no repetitive sequences detected" in lower_line
                or "only contains ambiguous bases" in lower_line
            ):
                break

            if line.startswith(("SW", "score", "There were")):
                continue

            try:
                parsed_row = parse_row(input_path, line)
            except ValueError as exc:
                errors.append(str(exc))
                continue

            if parsed_row is None:
                continue

            if parsed_row.consensus_triplet in consensus_keys_by_triplet:
                consensus_key = consensus_keys_by_triplet[parsed_row.consensus_triplet]
                parsed_row.feature["repeat_consensus"] = consensus_key

            features.append(parsed_row.feature)

    if errors:
        raise ValueError(format_parse_errors("RepeatMasker output", input_path, errors))

    return features, consensuses_by_key
