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
"""Unit testing of shared ``convert_to_genomio_json`` helpers."""

from contextlib import nullcontext as does_not_raise
from datetime import datetime, timezone
import hashlib
import logging
import os
from pathlib import Path
from typing import Callable, ContextManager

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import base


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


def test_consensus_sha256_key_normalises_fields() -> None:
    """Test ``convert_to_genomio_json.Consensus.sha256_key()`` normalises whitespace and sequence case."""
    consensus = convert_to_genomio_json.Consensus(
        name=" AluY ",
        repeat_class=" SINE ",
        repeat_type=" Alu ",
        seq="ac gt\n",
    )

    assert consensus.sha256_key() == _sha256_key("AluY", "SINE", "Alu", "ACGT")


@pytest.mark.parametrize(
    ("parser", "token", "field_name", "raw_line", "expectation"),
    [
        pytest.param(int, "42", "count", "42", does_not_raise(42), id="Valid integer"),
        pytest.param(
            int,
            "abc",
            "count",
            "abc",
            pytest.raises(ValueError, match=r"Invalid 'count' in input\.out: token='abc', line='abc'"),
            id="Invalid integer",
        ),
        pytest.param(float, "0.25", "score", "0.25", does_not_raise(0.25), id="Valid float"),
        pytest.param(
            float,
            "abc",
            "score",
            "abc",
            pytest.raises(ValueError, match=r"Invalid 'score' in input\.out: token='abc', line='abc'"),
            id="Invalid float",
        ),
    ],
)
def test_parse_token(
    parser: Callable[[str], int | float],
    token: str,
    field_name: str,
    raw_line: str,
    expectation: ContextManager,
) -> None:
    """Test ``base.parse_token()`` parses valid tokens and reports invalid token context.

    Args:
        parser: Parser callable to apply to the token.
        token: Raw token to parse.
        field_name: Field name used in error messages.
        raw_line: Original line used in error messages.
        expectation: Context manager for the expected result or exception.

    """
    with expectation as expected:
        assert base.parse_token(parser, token, field_name, raw_line, Path("input.out")) == expected


def test_format_parse_errors() -> None:
    """Test ``base.format_parse_errors()`` formats a counted bullet list."""
    assert base.format_parse_errors(
        "TRF output",
        Path("input.dat"),
        ["first error", "second error"],
    ) == ("Found 2 errors while parsing TRF output in input.dat:\n- first error\n- second error")


def test_file_last_modified_time_returns_utc_isoformat(tmp_path: Path) -> None:
    """Test ``base.file_last_modified_time()`` returns a UTC ISO timestamp.

    Args:
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = tmp_path / "input.out"
    input_path.write_text("content", encoding="utf-8")
    modified_time = datetime(2024, 1, 2, 3, 4, 5, tzinfo=timezone.utc).timestamp()

    os.utime(input_path, (modified_time, modified_time))

    assert base.file_last_modified_time(input_path) == "2024-01-02T03:04:05Z"


@pytest.mark.parametrize(
    ("seq_region_start", "seq_region_end", "repeat_start", "repeat_end", "expectation", "warning_pattern"),
    [
        pytest.param(1, 10, 2, 5, does_not_raise(enter_result=True), None, id="Valid coordinates"),
        pytest.param(
            0,
            10,
            1,
            5,
            pytest.raises(ValueError, match=r"Invalid seq_region coordinates"),
            None,
            id="Non-positive sequence region start",
        ),
        pytest.param(
            10,
            9,
            1,
            5,
            pytest.raises(ValueError, match=r"seq_region_end < seq_region_start"),
            None,
            id="Sequence region end before start",
        ),
        pytest.param(
            1,
            10,
            0,
            5,
            does_not_raise(enter_result=False),
            r"Invalid repeat coordinates",
            id="Non-positive repeat start",
        ),
        pytest.param(
            1,
            10,
            5,
            4,
            does_not_raise(enter_result=False),
            r"repeat_end < repeat_start",
            id="Repeat end before start",
        ),
    ],
)
def test_has_valid_parsed_coordinates(
    *,
    seq_region_start: int,
    seq_region_end: int,
    repeat_start: int,
    repeat_end: int,
    expectation: ContextManager,
    warning_pattern: str | None,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test ``base.has_valid_parsed_coordinates()`` correctly validates coordinates.

    Args:
        seq_region_start: Sequence region start coordinate.
        seq_region_end: Sequence region end coordinate.
        repeat_start: Repeat start coordinate.
        repeat_end: Repeat end coordinate.
        expectation: Context manager for the expected result or exception.
        warning_pattern: Optional substring expected to appear in the logged warning message.
        caplog: Pytest fixture for capturing log output.

    """
    with caplog.at_level(logging.WARNING), expectation as expected:
        assert (
            base.has_valid_parsed_coordinates(
                Path("input.out"),
                seq_region_start=seq_region_start,
                seq_region_end=seq_region_end,
                repeat_start=repeat_start,
                repeat_end=repeat_end,
                line="raw line",
            )
            == expected
        )

    if warning_pattern is None:
        assert not caplog.records
    else:
        assert any(warning_pattern in record.message for record in caplog.records)
