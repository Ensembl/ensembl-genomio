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

"""Unit testing of `ensembl.io.genomio.repeats.convert_to_genomio_json` module."""

from contextlib import nullcontext as does_not_raise
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
from typing import Callable, ContextManager
from unittest.mock import Mock, patch

import pytest
from pytest import param

from ensembl.io.genomio.features import convert_to_genomio_json


def _sha256_key(name: str, repeat_class: str, repeat_type: str, seq: str) -> str:
    """
    Computes the expected SHA-256 repeat consensus key.

    Args:
        name: Repeat name.
        repeat_class: Repeat class.
        repeat_type: Repeat type.
        seq: Consensus sequence.

    Returns:
        str:Expected SHA-256 digest.
    """
    payload = f"{name}\t{repeat_class}\t{repeat_type}\t{seq}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def test_consensus_sha256_key_normalises_fields() -> None:
    """
    Tests ``convert_to_genomio_json.Consensus.sha256_key()`` normalises whitespace and sequence case.
    """
    consensus = convert_to_genomio_json.Consensus(
        name=" AluY ",
        repeat_class=" SINE ",
        repeat_type=" Alu ",
        seq="ac gt\n",
    )

    assert consensus.sha256_key() == _sha256_key("AluY", "SINE", "Alu", "ACGT")


@pytest.mark.parametrize(
    ("repeat_class", "expected"),
    [
        ("LINE1", "Type I Transposons/LINE"),
        ("Satellite", "Satellite repeats"),
        ("fooRNA", "RNA repeats"),
        ("NotInMap", "Unknown"),
    ],
)
def test_map_repeatmasker_repeat_consensus_type(repeat_class: str, expected: str) -> None:
    """
    Tests the mapping of RepeatMasker repeat classes to GenomIO categories for known and unknown patterns.

    Args:
        repeat_class: Input repeat class string.
        expected: Expected mapped output.
    """
    assert convert_to_genomio_json._map_repeatmasker_repeat_consensus_type(repeat_class) == expected


@pytest.mark.parametrize(
    ("parser", "token", "field_name", "raw_line", "expectation"),
    [
        param(int, "42", "count", "42", does_not_raise(42), id="Valid integer"),
        param(
            int,
            "abc",
            "count",
            "abc",
            pytest.raises(ValueError, match=r"Invalid 'count' in input\.out: token='abc', line='abc'"),
            id="Invalid integer",
        ),
        param(float, "0.25", "score", "0.25", does_not_raise(0.25), id="Valid float"),
        param(
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
    """
    Tests ``convert_to_genomio_json._parse_token()`` parses valid tokens and reports invalid token context.

    Args:
        parser: Parser callable to apply to the token.
        token: Raw token to parse.
        field_name: Field name used in error messages.
        raw_line: Original line used in error messages.
        expectation: Context manager for the expected result or exception.
    """
    with expectation as expected:
        assert (
            convert_to_genomio_json._parse_token(parser, token, field_name, raw_line, Path("input.out"))
            == expected
        )


def test_format_parse_errors() -> None:
    """
    Tests ``convert_to_genomio_json._format_parse_errors()`` formats a counted bullet list.
    """
    assert convert_to_genomio_json._format_parse_errors(
        "TRF output",
        Path("input.dat"),
        ["first error", "second error"],
    ) == ("Found 2 errors while parsing TRF output in input.dat:\n" "- first error\n" "- second error")


def test_file_last_modified_time_returns_utc_isoformat(tmp_path: Path) -> None:
    """
    Tests ``convert_to_genomio_json._file_last_modified_time()`` returns a UTC ISO timestamp.

    Args:
        tmp_path: Temporary directory provided by pytest.
    """
    input_path = tmp_path / "input.out"
    input_path.write_text("content", encoding="utf-8")
    modified_time = datetime(2024, 1, 2, 3, 4, 5, tzinfo=timezone.utc).timestamp()

    os.utime(input_path, (modified_time, modified_time))

    assert convert_to_genomio_json._file_last_modified_time(input_path) == "2024-01-02T03:04:05Z"


@pytest.mark.parametrize(
    ("seq_region_start", "seq_region_end", "repeat_start", "repeat_end", "expectation", "warning_pattern"),
    [
        param(1, 10, 2, 5, does_not_raise(True), None, id="Valid coordinates"),
        param(
            0,
            10,
            1,
            5,
            pytest.raises(ValueError, match=r"Invalid seq_region coordinates"),
            None,
            id="Non-positive sequence region start",
        ),
        param(
            10,
            9,
            1,
            5,
            pytest.raises(ValueError, match=r"seq_region_end < seq_region_start"),
            None,
            id="Sequence region end before start",
        ),
        param(
            1,
            10,
            0,
            5,
            does_not_raise(False),
            r"Invalid repeat coordinates",
            id="Non-positive repeat start",
        ),
        param(
            1,
            10,
            5,
            4,
            does_not_raise(False),
            r"repeat_end < repeat_start",
            id="Repeat end before start",
        ),
    ],
)
def test_has_valid_parsed_coordinates(
    seq_region_start: int,
    seq_region_end: int,
    repeat_start: int,
    repeat_end: int,
    expectation: ContextManager,
    warning_pattern: str | None,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """
    Tests ``convert_to_genomio_json._has_valid_parsed_coordinates()`` correctly validates coordinates.

    Args:
        seq_region_start: Sequence region start coordinate.
        seq_region_end: Sequence region end coordinate.
        repeat_start: Repeat start coordinate.
        repeat_end: Repeat end coordinate.
        expectation: Context manager for the expected result or exception.
        warning_pattern: Optional substring expected to appear in the logged warning message.
        caplog: Pytest fixture for capturing log output.
    """
    with caplog.at_level(logging.WARNING):
        with expectation as expected:
            assert (
                convert_to_genomio_json._has_valid_parsed_coordinates(
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


@pytest.mark.parametrize(
    "entries",
    [
        [
            ("AluY", "SINE", "Alu", "Type I Transposons/SINE", "acgt"),
            ("Foo", "LINE", "Unknown", "Type I Transposons/LINE", "AATT"),
            ("Bar", "Unknown", "Unknown", "Unknown", "ccgg"),
        ]
    ],
)
def test_parse_repeatmasker_consensus_library(
    data_dir: Path,
    entries: list[tuple[str, str, str, str, str]],
) -> None:
    """
    Tests the ``convert_to_genomio_json._parse_repeatmasker_consensus_library()`` internal helper.

    Args:
        data_dir: Module's test data directory fixture.
        entries: Expected consensus name, class, raw repeat type, mapped repeat type, and sequence values.
    """
    fasta_path = data_dir / "repeatmodeler" / "classified_mixed.fa"
    expected_consensus_by_triplet = {
        (name, repeat_class, raw_repeat_type): convert_to_genomio_json.Consensus(
            name=name,
            repeat_class=repeat_class,
            repeat_type=expected_repeat_type,
            seq=seq,
        )
        for name, repeat_class, raw_repeat_type, expected_repeat_type, seq in entries
    }
    expected_consensuses = {
        consensus.sha256_key(): consensus for consensus in expected_consensus_by_triplet.values()
    }
    expected_keys_by_triplet = {
        triplet: consensus.sha256_key() for triplet, consensus in expected_consensus_by_triplet.items()
    }

    consensus_keys_by_triplet, consensuses_by_key = (
        convert_to_genomio_json._parse_repeatmasker_consensus_library(fasta_path)
    )

    assert consensus_keys_by_triplet == expected_keys_by_triplet
    assert consensuses_by_key == expected_consensuses


@pytest.mark.parametrize(
    ("filename", "expected_features"),
    [
        param(
            "repeatmasker_out/success_plus_with_star.out",
            [
                {
                    "seq_region": "chr1",
                    "seq_region_start": 10,
                    "seq_region_end": 20,
                    "seq_region_strand": "+",
                    "repeat_start": 1,
                    "repeat_end": 11,
                    "score": 100.0,
                    "attributes": {
                        "perc_div": 1.0,
                        "perc_del": 0.0,
                        "perc_ins": 0.0,
                        "repeatmasker_repeat_type": "Alu",
                    },
                }
            ],
            id="Forward strand entry with trailing *",
        ),
        param(
            "repeatmasker_out/success_c_strand.out",
            [
                {
                    "seq_region": "chr2",
                    "seq_region_start": 100,
                    "seq_region_end": 150,
                    "seq_region_strand": "-",
                    "repeat_start": 2,
                    "repeat_end": 40,
                    "score": 200.0,
                    "attributes": {
                        "perc_div": 2.5,
                        "perc_del": 1.5,
                        "perc_ins": 0.5,
                        "repeatmasker_repeat_type": "L1",
                    },
                }
            ],
            id="Complement strand entry",
        ),
        param(
            "repeatmasker_out/stop_no_repeats.out",
            [
                {
                    "seq_region": "chr3",
                    "seq_region_start": 7,
                    "seq_region_end": 9,
                    "seq_region_strand": "+",
                    "repeat_start": 3,
                    "repeat_end": 8,
                    "score": 123.0,
                    "attributes": {
                        "perc_div": 4.0,
                        "perc_del": 0.1,
                        "perc_ins": 0.2,
                        "repeatmasker_repeat_type": "Unknown",
                    },
                }
            ],
            id="Stops on no repeats message",
        ),
        param(
            "repeatmasker_out/stop_ambiguous_bases.out",
            [
                {
                    "seq_region": "chr4",
                    "seq_region_start": 3,
                    "seq_region_end": 6,
                    "seq_region_strand": "+",
                    "repeat_start": 2,
                    "repeat_end": 5,
                    "score": 321.0,
                    "attributes": {
                        "perc_div": 0.0,
                        "perc_del": 0.0,
                        "perc_ins": 0.0,
                        "repeatmasker_repeat_type": "Alu",
                    },
                }
            ],
            id="Stops on ambiguous bases in repeat coordinates",
        ),
        param(
            "repeatmasker_out/with_blank_lines.out",
            [
                {
                    "seq_region": "chr1",
                    "seq_region_start": 10,
                    "seq_region_end": 20,
                    "seq_region_strand": "+",
                    "repeat_start": 1,
                    "repeat_end": 11,
                    "score": 100.0,
                    "attributes": {
                        "perc_div": 1.0,
                        "perc_del": 0.0,
                        "perc_ins": 0.0,
                        "repeatmasker_repeat_type": "Alu",
                    },
                },
                {
                    "seq_region": "chr1",
                    "seq_region_start": 30,
                    "seq_region_end": 40,
                    "seq_region_strand": "+",
                    "repeat_start": 1,
                    "repeat_end": 11,
                    "score": 200.0,
                    "attributes": {
                        "perc_div": 2.0,
                        "perc_del": 0.0,
                        "perc_ins": 0.0,
                        "repeatmasker_repeat_type": "Unknown",
                    },
                },
            ],
            id="Skips blank lines",
        ),
    ],
)
def test_parse_repeatmasker_output_success(
    data_dir: Path,
    filename: str,
    expected_features: list[dict[str, object]],
) -> None:
    """
    Tests successful parsing of RepeatMasker `.out` files.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Input RepeatMasker `.out` filename.
        expected_features: Expected parsed repeat feature dictionaries.
    """
    out_path = data_dir / filename
    features, consensuses_by_key = convert_to_genomio_json.parse_repeatmasker_output(
        out_path,
        consensus_lib_path=None,
    )

    assert features == expected_features
    assert consensuses_by_key == {}


def test_parse_repeatmasker_output_adds_consensus_from_library(data_dir: Path) -> None:
    """
    Tests RepeatMasker parsing attaches consensus keys when a matching consensus library is provided.

    Args:
        data_dir: Module's test data directory fixture.
    """
    out_path = data_dir / "create_json" / "basic.out"
    consensus_lib_path = data_dir / "create_json" / "repeatmodeler_match.fa"
    expected_consensus = convert_to_genomio_json.Consensus(
        name="Foo",
        repeat_class="SINE",
        repeat_type="Type I Transposons/SINE",
        seq="acgt",
    )
    expected_consensus_key = expected_consensus.sha256_key()

    features, consensuses_by_key = convert_to_genomio_json.parse_repeatmasker_output(
        out_path,
        consensus_lib_path=consensus_lib_path,
    )

    assert features == [
        {
            "seq_region": "chr1",
            "seq_region_start": 10,
            "seq_region_end": 20,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 11,
            "repeat_consensus": expected_consensus_key,
            "score": 100.0,
            "attributes": {
                "perc_div": 1.0,
                "perc_del": 0.0,
                "perc_ins": 0.0,
                "repeatmasker_repeat_type": "Alu",
            },
        }
    ]
    assert consensuses_by_key == {expected_consensus_key: expected_consensus}


@pytest.mark.parametrize(
    ("filename", "error_pattern"),
    [
        param(
            "repeatmasker_out/error_13_columns.out",
            r"Expected 14 or 15 columns",
            id="Rejects lines with incorrect number of columns",
        ),
        param(
            "repeatmasker_out/error_invalid_score.out",
            r"Invalid 'score'",
            id="Rejects invalid score",
        ),
        param(
            "repeatmasker_out/error_invalid_perc_div.out",
            r"Invalid 'perc_div'",
            id="Rejects invalid perc_div",
        ),
        param(
            "repeatmasker_out/error_invalid_seq_region_start.out",
            r"Invalid 'seq_region_start'",
            id="Rejects invalid sequence region start",
        ),
        param(
            "repeatmasker_out/error_non_positive_seq_region_start.out",
            r"Invalid seq_region coordinates",
            id="Rejects non-positive sequence region start",
        ),
        param(
            "repeatmasker_out/error_seq_region_end_before_start.out",
            r"seq_region_end < seq_region_start",
            id="Rejects sequence region end before start",
        ),
        param(
            "repeatmasker_out/error_unexpected_strand.out",
            r"Unexpected strand token",
            id="Rejects unexpected strand token",
        ),
        param(
            "repeatmasker_out/error_empty_repeat_type.out",
            r"Malformed repeat_class/family",
            id="Rejects empty repeat type",
        ),
        param(
            "repeatmasker_out/error_empty_repeat_class.out",
            r"Malformed repeat_class/family",
            id="Rejects empty repeat class",
        ),
        param(
            "repeatmasker_out/error_plus_repeat_start_invalid.out",
            r"Invalid 'repeat_start'",
            id="Rejects invalid forward strand repeat start",
        ),
        param(
            "repeatmasker_out/error_c_repeat_start_invalid.out",
            r"Invalid 'repeat_start'",
            id="Rejects invalid reverse strand repeat start",
        ),
    ],
)
def test_parse_repeatmasker_output_errors(data_dir: Path, filename: str, error_pattern: str) -> None:
    """
    Tests that malformed RepeatMasker `.out` rows raise errors.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Name of RepeatMasker `.out` file containing invalid data.
        error_pattern: Regex expected in the raised error message.
    """
    out_path = data_dir / filename
    with pytest.raises(ValueError, match=error_pattern):
        convert_to_genomio_json.parse_repeatmasker_output(out_path, None)


@pytest.mark.parametrize(
    ("filename", "warning_pattern"),
    [
        param(
            "repeatmasker_out/error_non_positive_repeat_end.out",
            r"Invalid repeat coordinates",
            id="Skips non-positive repeat end with warning",
        ),
        param(
            "repeatmasker_out/error_repeat_end_before_start.out",
            r"repeat_end < repeat_start",
            id="Skips repeat end before start with warning",
        ),
    ],
)
def test_parse_repeatmasker_output_skips_invalid_repeat_coordinates(
    data_dir: Path,
    filename: str,
    warning_pattern: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """
    Tests that ``convert_to_genomio_json.parse_repeatmasker_output()`` logs a warning and skips invalid records.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Input RepeatMasker `.out` filename containing invalid coordinates.
        warning_pattern: Substring expected to appear in the logged warning message.
        caplog: Pytest fixture for capturing log output.
    """
    out_path = data_dir / filename

    with caplog.at_level("WARNING"):
        features, consensuses_by_key = convert_to_genomio_json.parse_repeatmasker_output(out_path, None)

    assert features == []
    assert consensuses_by_key == {}
    assert caplog.text
    assert any(warning_pattern in record.message for record in caplog.records)


@pytest.mark.parametrize(
    ("filename", "expected_features", "expected_consensuses"),
    [
        param(
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
        param(
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
    data_dir: Path,
    filename: str,
    expected_features: list[dict[str, object]],
    expected_consensuses: dict[str, convert_to_genomio_json.Consensus],
) -> None:
    """
    Tests successful parsing of TRF `.dat` examples.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Input TRF `.dat` filename.
        expected_features: Expected parsed repeat feature dictionaries.
        expected_consensuses: Expected parsed consensus records.
    """
    trf_path = data_dir / filename

    features, consensuses_by_key = convert_to_genomio_json.parse_trf_output(trf_path)

    assert features == expected_features
    assert consensuses_by_key == expected_consensuses


@pytest.mark.parametrize(
    ("parser_name", "filename", "expected_fragments"),
    [
        param(
            "repeatmasker",
            "repeatmasker_out/error_multiple_errors.out",
            [
                "Found 3 errors while parsing RepeatMasker output",
                "Invalid 'score'",
                "seq_region_end < seq_region_start",
                "Unexpected strand token",
            ],
            id="RepeatMasker multiple row errors",
        ),
        param(
            "trf",
            "trf/error_multiple_errors.dat",
            [
                "Found 3 errors while parsing TRF output",
                "Invalid 'start'",
                "seq_region_end < seq_region_start",
                "Invalid 'period_size'",
            ],
            id="TRF multiple row errors",
        ),
        param(
            "trf",
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
        param(
            "trf",
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
def test_parse_output_collates_all_errors(
    data_dir: Path,
    parser_name: str,
    filename: str,
    expected_fragments: list[str],
) -> None:
    """
    Tests that malformed parser rows are collated before raising.

    Args:
        data_dir: Module's test data directory fixture.
        parser_name: Parser to run against the input file.
        filename: Input filename containing invalid data.
        expected_fragments: Expected substrings in the raised error message.
    """
    input_path = data_dir / filename

    with pytest.raises(ValueError) as excinfo:
        if parser_name == "repeatmasker":
            convert_to_genomio_json.parse_repeatmasker_output(input_path, None)
        elif parser_name == "trf":
            convert_to_genomio_json.parse_trf_output(input_path)
        else:
            raise AssertionError(f"Unsupported parser name: {parser_name}")

    error_message = str(excinfo.value)
    for expected_fragment in expected_fragments:
        if expected_fragment.startswith("not:"):
            assert expected_fragment.removeprefix("not:") not in error_message
        else:
            assert expected_fragment in error_message


@patch("ensembl.io.genomio.features.convert_to_genomio_json.parse_repeatmasker_output")
def test_create_genomio_json_uses_repeatmasker_parser_output(
    mock_parse_repeatmasker_output: Mock,
    tmp_path: Path,
) -> None:
    """
    Tests JSON creation assembles mocked RepeatMasker parser output.

    Args:
        mock_parse_repeatmasker_output: Mock for ``convert_to_genomio_json.parse_repeatmasker_output()``.
        tmp_path: Temporary directory provided by pytest.
    """
    input = tmp_path / "input.out"
    input.write_text("parser input", encoding="utf-8")
    output = tmp_path / "out.json"
    consensus_lib = Path("consensus.fa")
    expected_features = [
        {
            "seq_region": "chr1",
            "seq_region_start": 10,
            "seq_region_end": 20,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 11,
            "repeat_consensus": "consensus-key",
            "score": 100.0,
            "attributes": {"repeatmasker_repeat_type": "Alu"},
        }
    ]
    expected_consensus = convert_to_genomio_json.Consensus(
        name="Foo",
        repeat_class="SINE",
        repeat_type="Type I Transposons/SINE",
        seq="acgt",
    )
    mock_parse_repeatmasker_output.return_value = (expected_features, {"consensus-key": expected_consensus})

    convert_to_genomio_json.create_genomio_json(
        input_path=input,
        repeatmasker_consensus_lib_path=consensus_lib,
        output_path=output,
        analysis_logic_name="repeatmask_customlib",
        analysis_display_label="Repeats: Custom library",
        analysis_description="desc",
        program="RepeatMasker",
        program_version="4.1.5",
        program_parameters="-nolow -gccalc -q",
        source_provider="Ensembl",
        is_primary=False,
    )

    mock_parse_repeatmasker_output.assert_called_once_with(input, consensus_lib)

    doc = json.loads(output.read_text(encoding="utf-8"))

    assert doc["analysis"]["logic_name"] == "repeatmask_customlib"
    assert doc["analysis"]["display_label"] == "Repeats: Custom library"
    assert doc["analysis"]["description"] == "desc"
    assert doc["analysis"]["program"] == "RepeatMasker"
    assert doc["analysis"]["program_version"] == "4.1.5"
    assert doc["analysis"]["program_parameters"] == "-nolow -gccalc -q"
    assert doc["analysis"]["run_date"].endswith("Z")

    assert doc["source"] == {"source_provider": "Ensembl", "is_primary": False}
    assert doc["repeat_features"] == expected_features
    assert doc["repeat_consensus"] == [
        {
            "repeat_consensus_key": "consensus-key",
            "repeat_name": "Foo",
            "repeat_class": "SINE",
            "repeat_type": "Type I Transposons/SINE",
            "repeat_consensus": "acgt",
        }
    ]


@patch("ensembl.io.genomio.features.convert_to_genomio_json.parse_trf_output")
def test_create_genomio_json_uses_trf_parser_output(
    mock_parse_trf_output: Mock,
    tmp_path: Path,
) -> None:
    """
    Tests JSON creation assembles mocked TRF parser output.

    Args:
        mock_parse_trf_output: Mock for ``convert_to_genomio_json.parse_trf_output()``.
        tmp_path: Temporary directory provided by pytest.
    """
    input = tmp_path / "input.out"
    input.write_text("parser input", encoding="utf-8")
    output = tmp_path / "out.json"
    expected_features = [
        {
            "seq_region": "chrA",
            "seq_region_start": 104,
            "seq_region_end": 109,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 2,
            "repeat_consensus": "trf-key",
            "score": 42.0,
            "attributes": {"period_size": 2},
        }
    ]
    expected_consensus = convert_to_genomio_json.Consensus(
        name="trf",
        repeat_class="trf",
        repeat_type="Tandem repeats",
        seq="AT",
    )
    mock_parse_trf_output.return_value = (expected_features, {"trf-key": expected_consensus})

    convert_to_genomio_json.create_genomio_json(
        input_path=input,
        repeatmasker_consensus_lib_path=None,
        output_path=output,
        analysis_logic_name="trf",
        analysis_display_label="Repeats: Custom library",
        analysis_description="desc",
        program="TRF",
        program_version="4.10.0",
        program_parameters="-h",
        source_provider="Ensembl",
        is_primary=False,
    )

    mock_parse_trf_output.assert_called_once_with(input)

    doc = json.loads(output.read_text(encoding="utf-8"))

    assert doc["analysis"]["logic_name"] == "trf"
    assert doc["analysis"]["program"] == "TRF"
    assert doc["analysis"]["program_version"] == "4.10.0"
    assert doc["analysis"]["program_parameters"] == "-h"
    assert doc["source"] == {"source_provider": "Ensembl", "is_primary": False}
    assert doc["repeat_features"] == expected_features
    assert doc["repeat_consensus"] == [
        {
            "repeat_consensus_key": "trf-key",
            "repeat_name": "trf",
            "repeat_class": "trf",
            "repeat_type": "Tandem repeats",
            "repeat_consensus": "AT",
        }
    ]


def test_create_genomio_json_rejects_unsupported_logic_name(data_dir: Path, tmp_path: Path) -> None:
    """
    Tests JSON creation rejects unsupported analysis logic names.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
    """
    with pytest.raises(ValueError, match=r"Unsupported analysis logic name"):
        convert_to_genomio_json.create_genomio_json(
            input_path=data_dir / "create_json" / "basic.out",
            repeatmasker_consensus_lib_path=None,
            output_path=tmp_path / "out.json",
            analysis_logic_name="unsupported",
            analysis_display_label="label",
            analysis_description="desc",
            program="RepeatMasker",
            program_version="4.1.5",
            program_parameters=None,
            source_provider="Ensembl",
            is_primary=False,
        )


@patch("ensembl.io.genomio.features.convert_to_genomio_json.parse_repeatmasker_output")
def test_create_genomio_json_omits_program_parameters_when_none(
    mock_parse_repeatmasker_output: Mock,
    tmp_path: Path,
) -> None:
    """
    Tests that program parameters are omitted when no value is provided.

    Args:
        mock_parse_repeatmasker_output: Mock for ``convert_to_genomio_json.parse_repeatmasker_output()``.
        tmp_path: Temporary directory provided by pytest.
    """
    input = tmp_path / "input.out"
    input.write_text("parser input", encoding="utf-8")
    output = tmp_path / "out.json"
    mock_parse_repeatmasker_output.return_value = ([], {})

    convert_to_genomio_json.create_genomio_json(
        input_path=input,
        repeatmasker_consensus_lib_path=None,
        output_path=output,
        analysis_logic_name="repeatmask_customlib",
        analysis_display_label="label",
        analysis_description="desc",
        program="RepeatMasker",
        program_version="4.1.5",
        program_parameters=None,
        source_provider="Ensembl",
        is_primary=True,
    )

    doc = json.loads(output.read_text(encoding="utf-8"))
    assert "program_parameters" not in doc["analysis"]
    assert doc["source"]["is_primary"] is True


@pytest.mark.parametrize(
    "use_repeatmodeler_lib",
    [
        param(False, id="Required arguments only"),
        param(True, id="All optional arguments overridden"),
    ],
)
def test_parse_args(data_dir: Path, tmp_path: Path, use_repeatmodeler_lib: bool) -> None:
    """
    Tests the ``convert_to_genomio_json.parse_args()`` function.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        use_repeatmodeler_lib: Whether to include optional repeatmodeler/program-parameters args.
    """
    input = data_dir / "create_json" / "basic.out"
    output = tmp_path / "out.json"

    argv = [
        "repeatmasker",
        "custom",
        "--input",
        str(input),
        "--output",
        str(output),
        "--program-version",
        "4.1.5",
    ]

    if use_repeatmodeler_lib:
        consensus_lib = data_dir / "create_json" / "repeatmodeler_match.fa"
        argv = [
            "repeatmasker",
            "repbase",
            "--input",
            str(input),
            "--consensus-lib",
            str(consensus_lib),
            "--output",
            str(output),
            "--program-version",
            "4.1.7",
            "--program-parameters",
            "-lib foo",
            "--source-provider",
            "Custom",
            "--is-primary",
        ]

    args = convert_to_genomio_json.parse_args(argv)

    assert args.__class__.__name__ == "Namespace"
    assert args.input == input
    assert args.output == output

    if use_repeatmodeler_lib:
        assert args.program_version == "4.1.7"
        assert args.analysis_logic_name == "repeatmask_rmlib"
        assert args.analysis_display_label == "Repeats: Repbase"
        assert (
            args.analysis_description
            == 'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, using the <a rel="external" href="http://www.girinst.org/repbase/">Repbase</a> library of repeat profiles.'
        )
        assert args.program == "RepeatMasker"
        assert args.source_provider == "Custom"
        assert args.is_primary is True
        assert args.consensus_lib == data_dir / "create_json" / "repeatmodeler_match.fa"
        assert args.program_parameters == "-lib foo"
    else:
        assert args.program_version == "4.1.5"
        assert args.analysis_logic_name == "repeatmask_customlib"
        assert args.analysis_display_label == "Repeats: Custom library"
        assert (
            args.analysis_description
            == 'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, using a custom library of <em>ab initio</em> repeat profiles for this species.'
        )
        assert args.program == "RepeatMasker"
        assert args.source_provider == "Ensembl"
        assert args.is_primary is False
        assert not hasattr(args, "consensus_lib")
        assert not hasattr(args, "program_parameters")


@patch("ensembl.io.genomio.features.convert_to_genomio_json.create_genomio_json")
@pytest.mark.parametrize(
    "use_repeatmodeler_lib",
    [
        param(False, id="Required arguments only, no consensus library"),
        param(True, id="All optional arguments overridden, consensus library used"),
    ],
)
def test_main_passes_expected_arguments(
    mock_create_genomio_json: Mock,
    data_dir: Path,
    tmp_path: Path,
    use_repeatmodeler_lib: bool,
) -> None:
    """
    Tests ``convert_to_genomio_json.main()`` passes parsed arguments correctly.

    Args:
        mock_create_genomio_json: Mock for ``convert_to_genomio_json.create_genomio_json()``.
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        use_repeatmodeler_lib: Whether to include optional repeatmodeler/program-parameters args.
    """
    input = data_dir / "create_json" / "basic.out"
    output = tmp_path / "out.json"

    argv = [
        "repeatmasker",
        "custom",
        "--input",
        str(input),
        "--output",
        str(output),
        "--program-version",
        "4.1.5",
    ]

    expected = {
        "input_path": input,
        "repeatmasker_consensus_lib_path": None,
        "output_path": output,
        "analysis_logic_name": "repeatmask_customlib",
        "analysis_display_label": "Repeats: Custom library",
        "analysis_description": 'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, using a custom library of <em>ab initio</em> repeat profiles for this species.',
        "program": "RepeatMasker",
        "program_version": "4.1.5",
        "program_parameters": None,
        "source_provider": "Ensembl",
        "is_primary": False,
    }

    if use_repeatmodeler_lib:
        consensus_lib = data_dir / "create_json" / "repeatmodeler_match.fa"
        argv = [
            "repeatmasker",
            "repbase",
            "--input",
            str(input),
            "--consensus-lib",
            str(consensus_lib),
            "--output",
            str(output),
            "--program-version",
            "4.1.7",
            "--program-parameters",
            "-lib foo",
            "--source-provider",
            "Custom",
            "--is-primary",
        ]
        expected = {
            "input_path": input,
            "repeatmasker_consensus_lib_path": consensus_lib,
            "output_path": output,
            "analysis_logic_name": "repeatmask_rmlib",
            "analysis_display_label": "Repeats: Repbase",
            "analysis_description": 'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, using the <a rel="external" href="http://www.girinst.org/repbase/">Repbase</a> library of repeat profiles.',
            "program": "RepeatMasker",
            "program_version": "4.1.7",
            "program_parameters": "-lib foo",
            "source_provider": "Custom",
            "is_primary": True,
        }

    mock_create_genomio_json.assert_called_once_with(**expected)


@patch("ensembl.io.genomio.features.convert_to_genomio_json.create_genomio_json")
def test_main_reraises_exceptions(mock_create_genomio_json: Mock, data_dir: Path, tmp_path: Path) -> None:
    """
    Tests the ``convert_to_genomio_json.main()`` function reraises exceptions.

    Args:
        mock_create_genomio_json: Mock for ``convert_to_genomio_json.create_genomio_json()``.
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
    """
    input = data_dir / "create_json" / "basic.out"
    output = tmp_path / "out.json"

    mock_create_genomio_json.side_effect = RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        convert_to_genomio_json.main(
            [
                "repeatmasker",
                "custom",
                "--input",
                str(input),
                "--output",
                str(output),
                "--program-version",
                "4.1.5",
            ]
        )
