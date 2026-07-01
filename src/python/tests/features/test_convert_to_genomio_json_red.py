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
"""Unit testing of Red GenomIO JSON conversion helpers."""

from pathlib import Path

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import red


def test_parse_row() -> None:
    """Test ``red.parse_row()`` parses one Red data row."""
    expected_consensus_key = red.RED_RPT_CONSENSUS.sha256_key()

    parsed_row = red.parse_row(Path("input.rpt"), "chr1 10 20")

    assert parsed_row.feature == {
        "seq_region": "chr1",
        "seq_region_start": 10,
        "seq_region_end": 20,
        "seq_region_strand": "+",
        "repeat_start": 1,
        "repeat_end": 11,
        "repeat_consensus": expected_consensus_key,
    }


@pytest.mark.parametrize(
    ("line", "error_pattern"),
    [
        pytest.param("chr1 10", r"Expected 3 columns", id="Too few columns"),
        pytest.param("chr1 10 20 extra", r"Expected 3 columns", id="Too many columns"),
        pytest.param("chr1 start 20", r"Invalid 'start'", id="Invalid start"),
        pytest.param("chr1 10 end", r"Invalid 'end'", id="Invalid end"),
        pytest.param("chr1 0 20", r"Invalid seq_region coordinates", id="Non-positive start"),
        pytest.param("chr1 10 0", r"Invalid seq_region coordinates", id="Non-positive end"),
        pytest.param("chr1 20 10", r"seq_region_end < seq_region_start", id="End before start"),
    ],
)
def test_parse_row_errors(line: str, error_pattern: str) -> None:
    """Test ``red.parse_row()`` rejects malformed Red data rows.

    Args:
        line: Raw Red data row.
        error_pattern: Regex expected in the raised error message.

    """
    with pytest.raises(ValueError, match=error_pattern):
        red.parse_row(Path("input.rpt"), line)


def test_parse_output_success(tmp_path: Path) -> None:
    """Test successful parsing of Red `.rpt` output.

    Args:
        tmp_path: Pytest temporary directory fixture.

    """
    rpt_path = tmp_path / "success.rpt"
    rpt_path.write_text("\nchr1 10 20\n\nchr2 30 30\n", encoding="utf-8")
    expected_consensus = red.RED_RPT_CONSENSUS
    expected_consensus_key = expected_consensus.sha256_key()

    features, consensuses_by_key = red.parse_output(rpt_path)

    assert features == [
        {
            "seq_region": "chr1",
            "seq_region_start": 10,
            "seq_region_end": 20,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 11,
            "repeat_consensus": expected_consensus_key,
        },
        {
            "seq_region": "chr2",
            "seq_region_start": 30,
            "seq_region_end": 30,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 1,
            "repeat_consensus": expected_consensus_key,
        },
    ]
    assert consensuses_by_key == {expected_consensus_key: expected_consensus}


def test_parse_output_collates_all_errors(tmp_path: Path) -> None:
    """Test that malformed Red parser rows are collated before raising.

    Args:
        tmp_path: Pytest temporary directory fixture.

    """
    rpt_path = tmp_path / "errors.rpt"
    rpt_path.write_text("chr1 start 20\nchr2 30 20\nchr3 40\n", encoding="utf-8")
    expected_fragments = [
        "Found 3 errors while parsing Red output",
        "Invalid 'start'",
        "seq_region_end < seq_region_start",
        "Expected 3 columns",
    ]

    with pytest.raises(
        ValueError,
        match=r"^Found \d+ errors while parsing Red output in .*:",
    ) as excinfo:
        red.parse_output(rpt_path)

    error_message = str(excinfo.value)
    for expected_fragment in expected_fragments:
        assert expected_fragment in error_message
