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
# pylint: disable=protected-access
"""Unit testing of TRF GenomIO JSON conversion helpers."""

from contextlib import nullcontext as does_not_raise
import hashlib
from pathlib import Path
from typing import ContextManager

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import trf


def _sha256_key(name: str, repeat_class: str, repeat_type: str, seq: str) -> str:
    """Compute the expected SHA-256 repeat consensus key.

    Args:
        name: Repeat name.
        repeat_class: Repeat class.
        repeat_type: Repeat type.
        seq: Consensus sequence.

    Returns:
        str:Expected SHA-256 digest.

    """
    payload = f"{name}\t{repeat_class}\t{repeat_type}\t{seq}".encode()
    return hashlib.sha256(payload).hexdigest()


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        pytest.param("Sequence: chrA:104-200", ("chrA", 104), id="Windowed sequence header"),
        pytest.param("Sequence: chrB", ("chrB", None), id="Plain sequence header"),
    ],
)
def test_parse_trf_sequence_header(line: str, expected: tuple[str, int | None]) -> None:
    """Test ``trf.parse_trf_sequence_header()`` extracts sequence ID and window start.

    Args:
        line: TRF sequence header line.
        expected: Expected sequence region and optional window start.

    """
    assert trf.parse_trf_sequence_header(line) == expected


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        pytest.param("Parameters: 2 7 7 80 10 50 500", "2 7 7 80 10 50 500", id="Parameters line"),
        pytest.param("Sequence: chrA", None, id="Different line type"),
    ],
)
def test_parse_trf_parameters(line: str, expected: str | None) -> None:
    """Test ``trf.parse_trf_parameters()`` extracts TRF parameters.

    Args:
        line: Raw TRF line.
        expected: Expected parameters string, or ``None`` for non-parameter lines.

    """
    assert trf.parse_trf_parameters(line) == expected


def test_missing_trf_sequence_error() -> None:
    """Test ``trf.missing_trf_sequence_error()`` includes skipped block context."""
    assert trf.missing_trf_sequence_error(Path("input.dat"), 3, 12) == (
        "TRF sequence header not found before data lines in input.dat: "
        "unprocessed_entries=3, first_data_line=12"
    )


@pytest.mark.parametrize(
    ("line", "expectation"),
    [
        pytest.param(
            "1 6 2 3.0 2 90.0 1.0 42 25 25 25 25 1.5 AT",
            does_not_raise(
                (
                    {
                        "seq_region": "chrA",
                        "seq_region_start": 104,
                        "seq_region_end": 109,
                        "seq_region_strand": "+",
                        "repeat_start": 1,
                        "repeat_end": 2,
                        "repeat_consensus": _sha256_key("trf", "trf", "Tandem repeats", "AT"),
                        "score": 42.0,
                        "attributes": {
                            "period_size": 2,
                            "copy_number": 3.0,
                            "consensus_size": 2,
                            "perc_match": 90.0,
                            "perc_indel": 1.0,
                            "entropy": 1.5,
                            "motif": "AT",
                            "a_pct": 25.0,
                            "c_pct": 25.0,
                            "g_pct": 25.0,
                            "t_pct": 25.0,
                            "trf_parameters": "2 7 7 80 10 50 500",
                        },
                    },
                    convert_to_genomio_json.Consensus(
                        name="trf",
                        repeat_class="trf",
                        repeat_type="Tandem repeats",
                        seq="AT",
                    ),
                )
            ),
            id="Valid data row",
        ),
        pytest.param(
            "1 6 2 3.0 2 90.0 1.0 42 25 25 25 25",
            pytest.raises(ValueError, match=r"Expected at least 13 columns"),
            id="Too few columns",
        ),
    ],
)
def test_parse_trf_data_row(
    line: str,
    expectation: ContextManager,
) -> None:
    """Test ``trf.parse_trf_data_row()`` parses or rejects one TRF data row.

    Args:
        line: Raw TRF data row.
        expectation: Context manager for the expected result or exception.

    """
    with expectation as expected:
        parsed_row = trf.parse_trf_data_row(
            Path("input.dat"),
            line,
            seq_region="chrA",
            window_start=104,
            trf_parameters="2 7 7 80 10 50 500",
        )

        expected_feature, expected_consensus = expected
        assert parsed_row.feature == expected_feature
        assert parsed_row.consensus == expected_consensus


@pytest.mark.parametrize(
    ("filename", "expected_features", "expected_consensuses"),
    [
        pytest.param(
            "trf/success_window_with_params.dat",
            [
                {
                    "seq_region": "chrA",
                    "seq_region_start": 104,
                    "seq_region_end": 109,
                    "seq_region_strand": "+",
                    "repeat_start": 1,
                    "repeat_end": 2,
                    "repeat_consensus": _sha256_key("trf", "trf", "Tandem repeats", "AT"),
                    "score": 42.0,
                    "attributes": {
                        "period_size": 2,
                        "copy_number": 3.0,
                        "consensus_size": 2,
                        "perc_match": 90.0,
                        "perc_indel": 1.0,
                        "entropy": 1.5,
                        "motif": "AT",
                        "a_pct": 25.0,
                        "c_pct": 25.0,
                        "g_pct": 25.0,
                        "t_pct": 25.0,
                        "trf_parameters": "2 7 7 80 10 50 500",
                    },
                }
            ],
            {
                _sha256_key("trf", "trf", "Tandem repeats", "AT"): convert_to_genomio_json.Consensus(
                    name="trf",
                    repeat_class="trf",
                    repeat_type="Tandem repeats",
                    seq="AT",
                )
            },
            id="Window coordinates with parameters",
        ),
        pytest.param(
            "trf/success_plain_no_params.dat",
            [
                {
                    "seq_region": "chrB",
                    "seq_region_start": 1,
                    "seq_region_end": 4,
                    "seq_region_strand": "+",
                    "repeat_start": 1,
                    "repeat_end": 4,
                    "repeat_consensus": _sha256_key("trf", "trf", "Tandem repeats", ""),
                    "score": 50.0,
                    "attributes": {
                        "period_size": 4,
                        "copy_number": 1.0,
                        "consensus_size": 4,
                        "perc_match": 100.0,
                        "perc_indel": 0.0,
                        "entropy": 2.0,
                        "motif": "",
                        "a_pct": 25.0,
                        "c_pct": 25.0,
                        "g_pct": 25.0,
                        "t_pct": 25.0,
                    },
                }
            ],
            {
                _sha256_key("trf", "trf", "Tandem repeats", ""): convert_to_genomio_json.Consensus(
                    name="trf",
                    repeat_class="trf",
                    repeat_type="Tandem repeats",
                    seq="",
                )
            },
            id="Plain coordinates without parameters or motif",
        ),
    ],
)
def test_parse_trf_output_success(
    convert_to_genomio_json_data_dir: Path,
    filename: str,
    expected_features: list[dict[str, object]],
    expected_consensuses: dict[str, convert_to_genomio_json.Consensus],
) -> None:
    """Test successful parsing of TRF `.dat` examples.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        filename: Input TRF `.dat` filename.
        expected_features: Expected parsed repeat feature dictionaries.
        expected_consensuses: Expected parsed consensus records.

    """
    trf_path = convert_to_genomio_json_data_dir / filename

    features, consensuses_by_key = convert_to_genomio_json.parse_trf_output(trf_path)

    assert features == expected_features
    assert consensuses_by_key == expected_consensuses


@pytest.mark.parametrize(
    ("filename", "expected_fragments"),
    [
        pytest.param(
            "trf/error_multiple_errors.dat",
            [
                "Found 3 errors while parsing TRF output",
                "Invalid 'start'",
                "seq_region_end < seq_region_start",
                "Invalid 'period_size'",
            ],
            id="TRF multiple row errors",
        ),
        pytest.param(
            "trf/error_missing_sequence_between_blocks.dat",
            [
                "Found 1 errors while parsing TRF output",
                "TRF sequence header not found before data lines",
                "unprocessed_entries=2",
                "first_data_line=7",
                "not:Invalid 'start'",
            ],
            id="TRF missing Sequence line between data blocks",
        ),
        pytest.param(
            "trf/error_missing_sequence_at_eof.dat",
            [
                "Found 1 errors while parsing TRF output",
                "TRF sequence header not found before data lines",
                "unprocessed_entries=2",
                "first_data_line=7",
                "not:Invalid 'start'",
            ],
            id="TRF missing Sequence line before EOF",
        ),
    ],
)
def test_parse_trf_output_collates_all_errors(
    convert_to_genomio_json_data_dir: Path,
    filename: str,
    expected_fragments: list[str],
) -> None:
    """Test that malformed TRF parser rows are collated before raising.

    Args:
        convert_to_genomio_json_data_dir: Test data directory for convert-to-GenomIO JSON fixtures.
        filename: Input filename containing invalid data.
        expected_fragments: Expected substrings in the raised error message.

    """
    input_path = convert_to_genomio_json_data_dir / filename

    with pytest.raises(
        ValueError,
        match=r"^Found \d+ errors while parsing TRF output in .*:",
    ) as excinfo:
        convert_to_genomio_json.parse_trf_output(input_path)

    error_message = str(excinfo.value)
    for expected_fragment in expected_fragments:
        if expected_fragment.startswith("not:"):
            assert expected_fragment.removeprefix("not:") not in error_message
        else:
            assert expected_fragment in error_message
