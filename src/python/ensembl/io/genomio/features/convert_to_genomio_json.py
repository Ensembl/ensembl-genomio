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
"""Constructs a GenomIO JSON document from output of feature identification tools."""

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import re
from typing import Callable, TypeVar

from Bio import SeqIO

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.archive import open_gz_file

__all__ = [
    "Consensus",
    "parse_repeatmasker_output",
    "parse_trf_output",
    "create_genomio_json",
]

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

T = TypeVar("T")

TRF_SEQUENCE_RE = re.compile(r"^Sequence:\s+(?P<seq_region>\S+?)(?::(?P<start>\d+)-\d+)?\s*$")

TRF_PARAMETERS_RE = re.compile(r"^Parameters:\s+(?P<params>.+)\s*$")


@dataclass(frozen=True)
class Consensus:
    """Repeat consensus record used for feature-to-consensus linking."""

    name: str
    repeat_class: str
    repeat_type: str
    seq: str

    def sha256_key(self) -> str:
        """
        Returns a normalized SHA256 digest for this consensus record.

        The digest is computed from the consensus name, repeat class, repeat type,
        and normalized sequence content.

        Returns:
            SHA256 hex digest for the consensus record.
        """
        norm_name = self.name.strip()
        norm_class = self.repeat_class.strip()
        norm_type = self.repeat_type.strip()
        norm_seq = "".join(self.seq.split()).upper()
        payload = f"{norm_name}\t{norm_class}\t{norm_type}\t{norm_seq}".encode("utf-8")
        return hashlib.sha256(payload).hexdigest()


def _parse_token(parser: Callable[[str], T], token: str, field_name: str, raw_line: str, path: Path) -> T:
    """
    Parses a field from tool output into the specified class type.

    Args:
        parser: A callable that takes a string and returns a value of type T.
        token: Raw token to parse.
        field_name: Field name used in error messages.
        raw_line: Original input line.
        path: Input file path.

    Returns:
        The parsed value.

    Raises:
        ValueError: If the token cannot be cast to the specified type.
    """
    try:
        return parser(token)
    except ValueError as exc:
        raise ValueError(f"Invalid {field_name!r} in {path}: token={token!r}, line={raw_line!r}") from exc


def _format_parse_errors(parser_name: str, input_path: Path, errors: list[str]) -> str:
    """Formats multiple parser errors into a single exception message."""
    return f"Found {len(errors)} errors while parsing {parser_name} in {input_path}:\n" + "\n".join(
        f"- {error}" for error in errors
    )


def _file_last_modified_time(file_path: Path) -> str:
    """Returns the last modified time of the given file."""
    return (
        datetime.fromtimestamp(
            file_path.stat().st_mtime,
            tz=timezone.utc,
        )
        .isoformat()
        .replace("+00:00", "Z")
    )


def _map_repeatmasker_repeat_consensus_type(repeat_class: str) -> str:
    """
    Maps a raw RepeatMasker repeat class to a GenomIO repeat category.

    Args:
        repeat_class: Raw repeat class string extracted from RepeatMasker output.

    Returns:
        Mapped repeat type category. Returns "Unknown" when no mapping matches.
    """
    for regex, mapped in REPEATMASKER_COMPILED_MAPPINGS:
        if regex.match(repeat_class):
            return mapped
    return "Unknown"


def _has_valid_parsed_coordinates(
    input_path: Path,
    *,
    seq_region_start: int,
    seq_region_end: int,
    repeat_start: int,
    repeat_end: int,
    line: str,
) -> bool:
    """
    Validates parsed coordinate values for a feature.

    Args:
        seq_region_start: Start coordinate on the sequence region.
        seq_region_end: End coordinate on the sequence region.
        repeat_start: Start coordinate on the repeat consensus.
        repeat_end: End coordinate on the repeat consensus.
        line: Original input line for error reporting.

    Returns:
        `True` if all coordinates are valid, `False` otherwise.

    Raises:
        ValueError: If any of the coordinate values are invalid (i.e. negative, zero, or end < start).
    """
    if seq_region_start < 1 or seq_region_end < 1:
        raise ValueError(
            f"Invalid seq_region coordinates in {input_path}: "
            f"start={seq_region_start}, end={seq_region_end}, line={line!r}"
        )
    if seq_region_end < seq_region_start:
        raise ValueError(
            f"seq_region_end < seq_region_start in {input_path}: "
            f"start={seq_region_start}, end={seq_region_end}, line={line!r}"
        )

    no_warnings = True
    if repeat_start < 1 or repeat_end < 1:
        logging.warning(
            f"Invalid repeat coordinates in {input_path}: "
            f"repeat_start={repeat_start}, repeat_end={repeat_end}, line={line!r}"
        )
        no_warnings = False
    if repeat_end < repeat_start:
        logging.warning(
            f"repeat_end < repeat_start in {input_path}: "
            f"repeat_start={repeat_start}, repeat_end={repeat_end}, line={line!r}"
        )
        no_warnings = False

    return no_warnings


def _parse_repeatmasker_consensus_library(
    consensus_lib_path: Path,
) -> tuple[dict[tuple[str, str, str], str], dict[str, Consensus]]:
    """
    Parses a RepeatMasker consensus library FASTA file into a dictionary of Consensus records.

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
                repeat_type=_map_repeatmasker_repeat_consensus_type(repeat_class),
                seq=str(record.seq),
            )
            consensus_key = consensus_obj.sha256_key()
            consensus_keys_by_triplet[(repeat_name, repeat_class, repeat_type)] = consensus_key
            consensuses_by_key[consensus_key] = consensus_obj
    return consensus_keys_by_triplet, consensuses_by_key


def parse_repeatmasker_output(
    input_path: Path, consensus_lib_path: Path | None
) -> tuple[list[dict], dict[str, Consensus]]:
    # pylint: disable=too-many-branches
    """
    Parse a RepeatMasker .out file into repeat feature dictionaries and consensus records.

    If a RepeatMasker consensus library FASTA file is provided, consensus sequences will be
    extracted from the library and associated with features based on their repeat names.
    Otherwise, consensus records will be created with empty sequences.

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
        consensus_keys_by_triplet, consensuses_by_key = _parse_repeatmasker_consensus_library(
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
                columns = line.split()
                if columns[-1] == "*":
                    columns.pop()

                if len(columns) < 14 or len(columns) > 15:
                    raise ValueError(
                        f"Expected 14 or 15 columns in {input_path}, got {len(columns)}: line={line!r}"
                    )

                score = _parse_token(float, columns[0], "score", line, input_path)
                perc_div = _parse_token(float, columns[1], "perc_div", line, input_path)
                perc_del = _parse_token(float, columns[2], "perc_del", line, input_path)
                perc_ins = _parse_token(float, columns[3], "perc_ins", line, input_path)

                seq_region = columns[4]
                seq_region_start = _parse_token(int, columns[5], "seq_region_start", line, input_path)
                seq_region_end = _parse_token(int, columns[6], "seq_region_end", line, input_path)

                strand_token = columns[8]
                if strand_token == "+":
                    seq_region_strand = "+"
                    repeat_start = _parse_token(int, columns[11], "repeat_start", line, input_path)
                    repeat_end = _parse_token(int, columns[12], "repeat_end", line, input_path)
                elif strand_token == "C":
                    seq_region_strand = "-"
                    repeat_end = _parse_token(int, columns[12], "repeat_end", line, input_path)
                    repeat_start = _parse_token(int, columns[13], "repeat_start", line, input_path)
                else:
                    raise ValueError(
                        f"Unexpected strand token in {input_path}: token={strand_token!r}, line={line!r}"
                    )

                repeat_name = columns[9]
                repeat_class_field = columns[10]

                if "/" in repeat_class_field:
                    repeat_class, repeat_type = repeat_class_field.split("/", 1)
                    if not repeat_class or not repeat_type:
                        raise ValueError(
                            f"Malformed repeat_class/family in {input_path}: "
                            f"value={repeat_class_field}, line={line!r}"
                        )
                else:
                    repeat_class = repeat_class_field
                    repeat_type = "Unknown"

                if not _has_valid_parsed_coordinates(
                    input_path,
                    seq_region_start=seq_region_start,
                    seq_region_end=seq_region_end,
                    repeat_start=repeat_start,
                    repeat_end=repeat_end,
                    line=line,
                ):
                    continue
            except ValueError as exc:
                errors.append(str(exc))
                continue

            feature = {
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
            }

            consensus_triplet = (repeat_name, repeat_class, repeat_type)
            if consensus_triplet in consensus_keys_by_triplet:
                consensus_key = consensus_keys_by_triplet[consensus_triplet]
                feature["repeat_consensus"] = consensus_key

            features.append(feature)

    if errors:
        raise ValueError(_format_parse_errors("RepeatMasker output", input_path, errors))

    return features, consensuses_by_key


def parse_trf_output(input_path: Path) -> tuple[list[dict], dict[str, Consensus]]:
    # pylint: disable=too-many-branches
    """
    Parses a TRF .dat file into repeat feature dictionaries and consensus records.

    Coordinates are converted from TRF's sequence-relative coordinates to genomic
    coordinates if the header provides a window like:
        Sequence: NC_135497.1:379500001-380000000

    Returns:
        A tuple containing:
            - repeat feature dictionaries
            - repeat consensus dictionary keyed by consensus SHA256 digest

    Raises:
        ValueError: If the TRF output contains malformed rows or invalid coordinate values.
    """
    features: list[dict] = []
    consensuses_by_key: dict[str, Consensus] = {}
    errors: list[str] = []

    seq_region: str | None = None
    window_start: int | None = None
    trf_parameters: str | None = None
    header_line = True
    in_data_block = False
    skip_data_block_without_sequence = False
    skipped_data_block_first_line: int | None = None
    skipped_data_block_entries = 0

    with open_gz_file(input_path) as fh:
        for line_number, raw_line in enumerate(fh, start=1):
            line = raw_line.strip()
            if not line:
                continue

            columns = line.split()
            if in_data_block and len(columns) < 13:
                if skip_data_block_without_sequence:
                    errors.append(
                        f"TRF sequence header not found before data lines in {input_path}: "
                        f"unprocessed_entries={skipped_data_block_entries}, "
                        f"first_data_line={skipped_data_block_first_line}"
                    )
                seq_region = None
                window_start = None
                trf_parameters = None
                in_data_block = False
                skip_data_block_without_sequence = False
                skipped_data_block_first_line = None
                skipped_data_block_entries = 0

            seq_match = TRF_SEQUENCE_RE.match(line)
            if seq_match is not None:
                header_line = False
                seq_region = seq_match.group("seq_region")
                if seq_match.group("start") is not None:
                    window_start = int(seq_match.group("start"))
                else:
                    window_start = None
                trf_parameters = None
                skip_data_block_without_sequence = False
                continue

            if header_line:
                continue

            params_match = TRF_PARAMETERS_RE.match(line)
            if params_match is not None:
                trf_parameters = params_match.group("params")
                continue

            if len(columns) < 13:
                continue

            if skip_data_block_without_sequence:
                skipped_data_block_entries += 1
                continue

            in_data_block = True

            try:
                if seq_region is None:
                    skip_data_block_without_sequence = True
                    skipped_data_block_first_line = line_number
                    skipped_data_block_entries = 1
                    continue

                start = _parse_token(int, columns[0], "start", line, input_path)
                end = _parse_token(int, columns[1], "end", line, input_path)
                period_size = _parse_token(int, columns[2], "period_size", line, input_path)
                copy_number = _parse_token(float, columns[3], "copy_number", line, input_path)
                consensus_size = _parse_token(int, columns[4], "consensus_size", line, input_path)
                perc_match = _parse_token(float, columns[5], "perc_match", line, input_path)
                perc_indel = _parse_token(float, columns[6], "perc_indel", line, input_path)
                score = _parse_token(float, columns[7], "score", line, input_path)
                a_pct = _parse_token(float, columns[8], "a_pct", line, input_path)
                c_pct = _parse_token(float, columns[9], "c_pct", line, input_path)
                g_pct = _parse_token(float, columns[10], "g_pct", line, input_path)
                t_pct = _parse_token(float, columns[11], "t_pct", line, input_path)
                entropy = _parse_token(float, columns[12], "entropy", line, input_path)
                motif = columns[13] if len(columns) >= 14 else ""

                if window_start is not None:
                    seq_region_start = window_start + start - 1
                    seq_region_end = window_start + end - 1
                else:
                    seq_region_start = start
                    seq_region_end = end

                _has_valid_parsed_coordinates(
                    input_path,
                    seq_region_start=seq_region_start,
                    seq_region_end=seq_region_end,
                    repeat_start=1,
                    repeat_end=period_size,
                    line=line,
                )
            except ValueError as exc:
                errors.append(str(exc))
                continue

            repeat_consensus = Consensus(
                name="trf",
                repeat_class="trf",
                repeat_type="Tandem repeats",
                seq="".join(motif.split()).upper(),
            )
            repeat_consensus_key = repeat_consensus.sha256_key()
            consensuses_by_key[repeat_consensus_key] = repeat_consensus

            attributes: dict[str, object] = {
                "period_size": period_size,
                "copy_number": copy_number,
                "consensus_size": consensus_size,
                "perc_match": perc_match,
                "perc_indel": perc_indel,
                "entropy": entropy,
                "motif": motif,
                "a_pct": a_pct,
                "c_pct": c_pct,
                "g_pct": g_pct,
                "t_pct": t_pct,
            }

            if trf_parameters is not None:
                attributes["trf_parameters"] = trf_parameters

            feature = {
                "seq_region": seq_region,
                "seq_region_start": seq_region_start,
                "seq_region_end": seq_region_end,
                "seq_region_strand": "+",
                "repeat_start": 1,
                "repeat_end": period_size,
                "repeat_consensus": repeat_consensus.sha256_key(),
                "score": score,
                "attributes": attributes,
            }

            features.append(feature)

    if skip_data_block_without_sequence:
        errors.append(
            f"TRF sequence header not found before data lines in {input_path}: "
            f"unprocessed_entries={skipped_data_block_entries}, "
            f"first_data_line={skipped_data_block_first_line}"
        )

    if errors:
        raise ValueError(_format_parse_errors("TRF output", input_path, errors))

    return features, consensuses_by_key


def create_genomio_json(
    input_path: Path,
    output_path: Path,
    *,
    analysis_logic_name: str,
    analysis_display_label: str,
    analysis_description: str,
    program: str,
    program_version: str,
    source_provider: str,
    is_primary: bool,
    repeatmasker_consensus_lib_path: Path | None = None,
    program_parameters: str | None = None,
) -> None:
    """
    Creates a GenomIO JSON document from feature identification tool output.

    Args:
        input: Path to the input file containing repeat masking results (e.g. RepeatMasker .out file).
        output: Path to the output JSON file to be created.
        analysis_logic_name: Logic name for the analysis (e.g. "repeatmask_customlib").
        analysis_display_label: Display label for the analysis.
        analysis_description: Description of the analysis (HTML allowed).
        program: Name of the program used to identify the features.
        program_version: Version of the program used to identify the features.
        source_provider: Name of the source provider for the features (e.g. "Ensembl").
        is_primary: Whether these features are primary or secondary annotations.
        repeatmasker_consensus_lib: Optional path to a FASTA file containing consensus sequences for the
            RepeatMasker library used.
        program_parameters: Optional parameters supplied to the program used for feature identification.

    Raises:
        ValueError: If an unsupported analysis logic name is provided.
    """
    features: list[dict[str, object]] = []
    consensuses_by_key: dict[str, Consensus] = {}

    if analysis_logic_name.startswith("repeatmask_"):
        features, consensuses_by_key = parse_repeatmasker_output(input_path, repeatmasker_consensus_lib_path)
    elif analysis_logic_name == "trf":
        features, consensuses_by_key = parse_trf_output(input_path)
    else:
        raise ValueError(f"Unsupported analysis logic name: {analysis_logic_name}")

    # Create dictionary for analysis metadata
    analysis: dict[str, str] = {
        "run_date": _file_last_modified_time(input_path),
        "logic_name": analysis_logic_name,
        "display_label": analysis_display_label,
        "description": analysis_description,
        "program": program,
        "program_version": program_version,
    }
    if program_parameters is not None:
        analysis["program_parameters"] = program_parameters

    # Construct JSON document
    json_doc: dict[str, object] = {
        "analysis": analysis,
        "source": {
            "source_provider": source_provider,
            "is_primary": bool(is_primary),
        },
        "repeat_features": features,
    }

    # Add consensus sequences if available
    if consensuses_by_key:
        repeat_consensuses: list[dict[str, str]] = []
        for consensus_key, consensus in consensuses_by_key.items():
            repeat_consensuses.append(
                {
                    "repeat_consensus_key": consensus_key,
                    "repeat_name": consensus.name,
                    "repeat_class": consensus.repeat_class,
                    "repeat_type": consensus.repeat_type,
                    "repeat_consensus": consensus.seq,
                }
            )
        json_doc["repeat_consensus"] = repeat_consensuses

    # Write JSON document to output file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(json_doc, indent=2) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments for JSON conversion.

    Args:
        argv: Optional list of command-line arguments. If `None`, arguments are
            taken from ``sys.argv``.

    Returns:
        Parsed command-line arguments.
    """

    def _add_common_arguments(subparser: ArgumentParser) -> None:

        subparser.add_argument_src_path(
            "--input",
            metavar="IN",
            required=True,
            help="Input file to be converted.",
        )
        subparser.add_argument_dst_path("--output", metavar="JSON", required=True, help="Output JSON path.")
        subparser.add_argument(
            "--program-version",
            required=True,
            help="Version of the program used to identify the features.",
        )
        subparser.add_argument(
            "--program-parameters",
            default=argparse.SUPPRESS,
            help="Parameters supplied to the program used for feature identification.",
        )
        subparser.add_argument(
            "--source-provider",
            default="Ensembl",
            help="Source provider for the features (default: Ensembl).",
        )
        subparser.add_argument(
            "--is-primary",
            action="store_true",
            help="Whether the source provider is the primary source for these features.",
        )
        subparser.add_log_arguments()

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)

    subparsers = parser.add_subparsers(dest="tool", required=True)

    # TRF subparser
    trf_parser = subparsers.add_parser("trf", help="Convert TRF output to GenomIO JSON.")
    _add_common_arguments(trf_parser)
    trf_parser.set_defaults(
        analysis_logic_name="trf",
        analysis_display_label="Tandem repeats (TRF)",
        analysis_description=(
            '<a rel="external" href="https://tandem.bu.edu/trf/trf.html">Tandem Repeats Finder</a> '
            "locates adjacent copies of a pattern of nucleotides."
        ),
        program="trf",
        repeatmasker_consensus_lib_path=None,
    )

    # RepeatMasker
    repeatmasker_parser = subparsers.add_parser(
        "repeatmasker", help="Convert RepeatMasker output to GenomIO JSON."
    )
    repeatmasker_parser.set_defaults(program="RepeatMasker")
    repeatmasker_subparsers = repeatmasker_parser.add_subparsers(dest="repeatmasker_mode", required=True)

    def _add_repeatmasker_common_args(subparser: ArgumentParser) -> None:
        _add_common_arguments(subparser)
        subparser.add_argument_src_path(
            "--consensus-lib",
            metavar="RM_LIB",
            default=argparse.SUPPRESS,
            help="FASTA file containing consensus sequences for the RepeatMasker library used.",
        )

    # RepeatMasker / custom
    custom_parser = repeatmasker_subparsers.add_parser(
        "custom",
        help="Convert RepeatMasker output generated using a custom library.",
    )
    _add_repeatmasker_common_args(custom_parser)
    custom_parser.set_defaults(
        analysis_logic_name="repeatmask_customlib",
        analysis_display_label="Repeats: Custom library",
        analysis_description=(
            'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, '
            "using a custom library of <em>ab initio</em> repeat profiles for this species."
        ),
    )

    # RepeatMasker / repbase
    repbase_parser = repeatmasker_subparsers.add_parser(
        "repbase",
        help="Convert RepeatMasker output geneterated using Repbase.",
    )
    _add_repeatmasker_common_args(repbase_parser)
    repbase_parser.set_defaults(
        analysis_logic_name="repeatmask_rmlib",
        analysis_display_label="Repeats: Repbase",
        analysis_description=(
            'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, '
            'using the <a rel="external" href="http://www.girinst.org/repbase/">Repbase</a> library of '
            "repeat profiles."
        ),
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """
    Run the JSON conversion command-line entry point.

    Args:
        argv: Optional list of command-line arguments. If `None`, arguments are taken from ``sys.argv``.
    """
    args = parse_args(argv)
    try:
        create_genomio_json(
            input_path=args.input,
            output_path=args.output,
            analysis_logic_name=args.analysis_logic_name,
            analysis_display_label=args.analysis_display_label,
            analysis_description=args.analysis_description,
            program=args.program,
            program_version=args.program_version,
            source_provider=args.source_provider,
            is_primary=args.is_primary,
            repeatmasker_consensus_lib_path=getattr(args, "consensus_lib", None),
            program_parameters=getattr(args, "program_parameters", None),
        )
    except Exception:
        logging.exception(f"Error processing file {args.input}")
        raise
