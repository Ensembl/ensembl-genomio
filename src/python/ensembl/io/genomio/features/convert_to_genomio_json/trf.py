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
"""Parse TRF output into GenomIO repeat feature records."""

from dataclasses import dataclass
from pathlib import Path
import re

from ensembl.io.genomio.features.convert_to_genomio_json.base import (
    Consensus,
    format_parse_errors,
    has_valid_parsed_coordinates,
    parse_token,
)
from ensembl.utils.archive import open_gz_file

MIN_TRF_COLUMNS = 13

TRF_SEQUENCE_RE = re.compile(r"^Sequence:\s+(?P<seq_region>\S+?)(?::(?P<start>\d+)-\d+)?\s*$")

TRF_PARAMETERS_RE = re.compile(r"^Parameters:\s+(?P<params>.+)\s*$")

__all__ = [
    "TRFParsedRow",
    "missing_trf_sequence_error",
    "parse_trf_data_row",
    "parse_trf_output",
    "parse_trf_parameters",
    "parse_trf_sequence_header",
]


@dataclass(frozen=True)
class TRFParsedRow:
    """Parsed TRF row and its repeat consensus record."""

    feature: dict[str, object]
    consensus: Consensus


def parse_trf_sequence_header(line: str) -> tuple[str, int | None] | None:
    """Parse a TRF sequence header line.

    Args:
        line: Raw TRF line without surrounding whitespace.

    Returns:
        Sequence region and optional window start, or ``None`` if the line is not a sequence header.

    """
    seq_match = TRF_SEQUENCE_RE.match(line)
    if seq_match is None:
        return None

    seq_region = seq_match.group("seq_region")
    window_start = int(seq_match.group("start")) if seq_match.group("start") is not None else None
    return seq_region, window_start


def parse_trf_parameters(line: str) -> str | None:
    """Parse a TRF parameters line, returning ``None`` if the line is not a parameters line.

    Args:
        line: Raw TRF line without surrounding whitespace.

    Returns:
        TRF parameters string, or ``None`` if the line is not a parameters line.

    """
    params_match = TRF_PARAMETERS_RE.match(line)
    return params_match.group("params") if params_match is not None else None


def missing_trf_sequence_error(
    input_path: Path,
    skipped_data_block_entries: int,
    skipped_data_block_first_line: int | None,
) -> str:
    """Format an error for TRF data rows that were not preceded by a Sequence header.

    Args:
        input_path: Input path of processed file.
        skipped_data_block_entries: Number of skipped data block entries.
        skipped_data_block_first_line: Line number of the first skipped data block line.

    Returns:
        Formatted error message describing the missing Sequence header and skipped data block.

    """
    return (
        f"TRF sequence header not found before data lines in {input_path}: "
        f"unprocessed_entries={skipped_data_block_entries}, "
        f"first_data_line={skipped_data_block_first_line}"
    )


def parse_trf_data_row(
    input_path: Path,
    line: str,
    *,
    seq_region: str,
    window_start: int | None,
    trf_parameters: str | None,
) -> TRFParsedRow:
    """Parse a single TRF data row.

    Args:
        input_path: Input TRF output path used for error messages.
        line: Raw TRF data row without surrounding whitespace.
        seq_region: Current sequence region from the preceding Sequence header.
        window_start: Optional genomic window start from the preceding Sequence header.
        trf_parameters: Optional current TRF parameters string.

    Returns:
        Parsed row containing feature data and consensus.

    Raises:
        ValueError: If the row is malformed or contains invalid coordinates.

    """
    columns = line.split()
    if len(columns) < MIN_TRF_COLUMNS:
        raise ValueError(
            f"Expected at least {MIN_TRF_COLUMNS} columns in {input_path}, got {len(columns)}: line={line!r}"
        )

    start = parse_token(int, columns[0], "start", line, input_path)
    end = parse_token(int, columns[1], "end", line, input_path)
    period_size = parse_token(int, columns[2], "period_size", line, input_path)
    copy_number = parse_token(float, columns[3], "copy_number", line, input_path)
    consensus_size = parse_token(int, columns[4], "consensus_size", line, input_path)
    perc_match = parse_token(float, columns[5], "perc_match", line, input_path)
    perc_indel = parse_token(float, columns[6], "perc_indel", line, input_path)
    score = parse_token(float, columns[7], "score", line, input_path)
    a_pct = parse_token(float, columns[8], "a_pct", line, input_path)
    c_pct = parse_token(float, columns[9], "c_pct", line, input_path)
    g_pct = parse_token(float, columns[10], "g_pct", line, input_path)
    t_pct = parse_token(float, columns[11], "t_pct", line, input_path)
    entropy = parse_token(float, columns[12], "entropy", line, input_path)
    motif = columns[13] if len(columns) >= 14 else ""  # noqa: PLR2004

    if window_start is not None:
        seq_region_start = window_start + start - 1
        seq_region_end = window_start + end - 1
    else:
        seq_region_start = start
        seq_region_end = end

    has_valid_parsed_coordinates(
        input_path,
        seq_region_start=seq_region_start,
        seq_region_end=seq_region_end,
        repeat_start=1,
        repeat_end=period_size,
        line=line,
    )

    repeat_consensus = Consensus(
        name="trf",
        repeat_class="trf",
        repeat_type="Tandem repeats",
        seq="".join(motif.split()).upper(),
    )

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

    return TRFParsedRow(
        feature={
            "seq_region": seq_region,
            "seq_region_start": seq_region_start,
            "seq_region_end": seq_region_end,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": period_size,
            "repeat_consensus": repeat_consensus.sha256_key(),
            "score": score,
            "attributes": attributes,
        },
        consensus=repeat_consensus,
    )


def parse_trf_output(input_path: Path) -> tuple[list[dict], dict[str, Consensus]]:  # noqa: PLR0912, PLR0915
    """Parse a TRF .dat file into repeat feature dictionaries and consensus records.

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
            if in_data_block and len(columns) < MIN_TRF_COLUMNS:
                if skip_data_block_without_sequence:
                    errors.append(
                        missing_trf_sequence_error(
                            input_path,
                            skipped_data_block_entries,
                            skipped_data_block_first_line,
                        )
                    )
                seq_region = None
                window_start = None
                trf_parameters = None
                in_data_block = False
                skip_data_block_without_sequence = False
                skipped_data_block_first_line = None
                skipped_data_block_entries = 0

            sequence_header = parse_trf_sequence_header(line)
            if sequence_header is not None:
                header_line = False
                seq_region, window_start = sequence_header
                trf_parameters = None
                skip_data_block_without_sequence = False
                continue

            if header_line:
                continue

            parsed_parameters = parse_trf_parameters(line)
            if parsed_parameters is not None:
                trf_parameters = parsed_parameters
                continue

            if len(columns) < MIN_TRF_COLUMNS:
                continue

            if skip_data_block_without_sequence:
                skipped_data_block_entries += 1
                continue

            in_data_block = True

            if seq_region is None:
                skip_data_block_without_sequence = True
                skipped_data_block_first_line = line_number
                skipped_data_block_entries = 1
                continue

            try:
                parsed_row = parse_trf_data_row(
                    input_path,
                    line,
                    seq_region=seq_region,
                    window_start=window_start,
                    trf_parameters=trf_parameters,
                )
            except ValueError as exc:
                errors.append(str(exc))
                continue

            consensuses_by_key[parsed_row.consensus.sha256_key()] = parsed_row.consensus
            features.append(parsed_row.feature)

    if skip_data_block_without_sequence:
        errors.append(
            missing_trf_sequence_error(
                input_path,
                skipped_data_block_entries,
                skipped_data_block_first_line,
            )
        )

    if errors:
        raise ValueError(format_parse_errors("TRF output", input_path, errors))

    return features, consensuses_by_key


_parse_trf_sequence_header = parse_trf_sequence_header
_parse_trf_parameters = parse_trf_parameters
_missing_trf_sequence_error = missing_trf_sequence_error
_parse_trf_data_row = parse_trf_data_row
