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
"""Unit testing of RepeatMasker GenomIO JSON conversion helpers."""

from pathlib import Path

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import repeatmasker


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
    """Test the mapping of RepeatMasker repeat classes to GenomIO categories for known and unknown patterns.

    Args:
        repeat_class: Input repeat class string.
        expected: Expected mapped output.

    """
    assert repeatmasker.map_repeatmasker_repeat_consensus_type(repeat_class) == expected


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
    convert_to_genomio_json_data_dir: Path,
    entries: list[tuple[str, str, str, str, str]],
) -> None:
    """Test the ``repeatmasker.parse_repeatmasker_consensus_library()`` internal helper.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        entries: Expected consensus name, class, raw repeat type, mapped repeat type, and sequence values.

    """
    fasta_path = convert_to_genomio_json_data_dir / "repeatmodeler" / "classified_mixed.fa"
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

    consensus_keys_by_triplet, consensuses_by_key = repeatmasker.parse_repeatmasker_consensus_library(
        fasta_path
    )

    assert consensus_keys_by_triplet == expected_keys_by_triplet
    assert consensuses_by_key == expected_consensuses


@pytest.mark.parametrize(
    ("repeat_class_field", "expected"),
    [
        pytest.param("SINE/Alu", ("SINE", "Alu"), id="Class and type"),
        pytest.param("Simple_repeat", ("Simple_repeat", "Unknown"), id="Class only"),
    ],
)
def test_parse_repeat_class_field(
    repeat_class_field: str,
    expected: tuple[str, str],
) -> None:
    """Test ``repeatmasker.parse_repeat_class_field()`` splits class fields.

    Args:
        repeat_class_field: Raw RepeatMasker repeat class/family field.
        expected: Expected repeat class and raw repeat type.

    """
    output = repeatmasker.parse_repeat_class_field(Path("input.out"), repeat_class_field, "raw line")
    assert output == expected


@pytest.mark.parametrize(
    "repeat_class_field",
    [
        pytest.param("SINE/", id="Missing type"),
        pytest.param("/Alu", id="Missing class"),
    ],
)
def test_parse_repeat_class_field_errors(repeat_class_field: str) -> None:
    """Test ``repeatmasker.parse_repeat_class_field()`` rejects malformed fields.

    Args:
        repeat_class_field: Raw RepeatMasker repeat class/family field.

    """
    with pytest.raises(ValueError, match=r"Malformed repeat_class/family"):
        repeatmasker.parse_repeat_class_field(
            Path("input.out"),
            repeat_class_field,
            "raw line",
        )


@pytest.mark.parametrize(
    ("columns", "expected"),
    [
        pytest.param(
            ["100", "1.0", "0.0", "0.0", "chr1", "10", "20", "(0)", "+", "Foo", "SINE/Alu", "1", "11", "(0)"],
            ("+", 1, 11),
            id="Forward strand",
        ),
        pytest.param(
            [
                "200",
                "2.5",
                "1.5",
                "0.5",
                "chr2",
                "100",
                "150",
                "(0)",
                "C",
                "Foo",
                "LINE/L1",
                "(0)",
                "40",
                "2",
            ],
            ("-", 2, 40),
            id="Complement strand",
        ),
    ],
)
def test_parse_strand_coordinates(
    columns: list[str],
    expected: tuple[str, int, int],
) -> None:
    """Test ``repeatmasker.parse_strand_coordinates()`` parses strand coordinates.

    Args:
        columns: Split RepeatMasker data row columns.
        expected: Expected sequence-region strand, repeat start, and repeat end.

    """
    assert (
        repeatmasker.parse_strand_coordinates(
            Path("input.out"),
            columns,
            "raw line",
        )
        == expected
    )


def test_parse_row() -> None:
    """Test ``repeatmasker.parse_row()`` parses one RepeatMasker data row."""
    line = "100 1.0 0.0 0.0 chr1 10 20 (0) + Foo SINE/Alu 1 11 (0) *"

    parsed_row = repeatmasker.parse_row(Path("input.out"), line)

    assert parsed_row is not None
    assert parsed_row.feature == {
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
    assert parsed_row.consensus_triplet == ("Foo", "SINE", "Alu")


@pytest.mark.parametrize(
    ("filename", "expected_features"),
    [
        pytest.param(
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
        pytest.param(
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
        pytest.param(
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
        pytest.param(
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
        pytest.param(
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
def test_parse_output_success(
    convert_to_genomio_json_data_dir: Path,
    filename: str,
    expected_features: list[dict[str, object]],
) -> None:
    """Test successful parsing of RepeatMasker `.out` files.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        filename: Input RepeatMasker `.out` filename.
        expected_features: Expected parsed repeat feature dictionaries.

    """
    out_path = convert_to_genomio_json_data_dir / filename
    features, consensuses_by_key = repeatmasker.parse_output(
        out_path,
        consensus_lib_path=None,
    )

    assert features == expected_features
    assert not consensuses_by_key


def test_parse_output_adds_consensus_from_library(
    convert_to_genomio_json_data_dir: Path,
) -> None:
    """Test RepeatMasker parsing attaches consensus keys when a matching consensus library is provided.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.

    """
    out_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    consensus_lib_path = convert_to_genomio_json_data_dir / "create_json" / "repeatmodeler_match.fa"
    expected_consensus = convert_to_genomio_json.Consensus(
        name="Foo",
        repeat_class="SINE",
        repeat_type="Type I Transposons/SINE",
        seq="acgt",
    )
    expected_consensus_key = expected_consensus.sha256_key()

    features, consensuses_by_key = repeatmasker.parse_output(
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
        pytest.param(
            "repeatmasker_out/error_13_columns.out",
            r"Expected 14 or 15 columns",
            id="Rejects lines with incorrect number of columns",
        ),
        pytest.param(
            "repeatmasker_out/error_invalid_score.out",
            r"Invalid 'score'",
            id="Rejects invalid score",
        ),
        pytest.param(
            "repeatmasker_out/error_invalid_perc_div.out",
            r"Invalid 'perc_div'",
            id="Rejects invalid perc_div",
        ),
        pytest.param(
            "repeatmasker_out/error_invalid_seq_region_start.out",
            r"Invalid 'seq_region_start'",
            id="Rejects invalid sequence region start",
        ),
        pytest.param(
            "repeatmasker_out/error_non_positive_seq_region_start.out",
            r"Invalid seq_region coordinates",
            id="Rejects non-positive sequence region start",
        ),
        pytest.param(
            "repeatmasker_out/error_seq_region_end_before_start.out",
            r"seq_region_end < seq_region_start",
            id="Rejects sequence region end before start",
        ),
        pytest.param(
            "repeatmasker_out/error_unexpected_strand.out",
            r"Unexpected strand token",
            id="Rejects unexpected strand token",
        ),
        pytest.param(
            "repeatmasker_out/error_empty_repeat_type.out",
            r"Malformed repeat_class/family",
            id="Rejects empty repeat type",
        ),
        pytest.param(
            "repeatmasker_out/error_empty_repeat_class.out",
            r"Malformed repeat_class/family",
            id="Rejects empty repeat class",
        ),
        pytest.param(
            "repeatmasker_out/error_plus_repeat_start_invalid.out",
            r"Invalid 'repeat_start'",
            id="Rejects invalid forward strand repeat start",
        ),
        pytest.param(
            "repeatmasker_out/error_c_repeat_start_invalid.out",
            r"Invalid 'repeat_start'",
            id="Rejects invalid reverse strand repeat start",
        ),
    ],
)
def test_parse_output_errors(
    convert_to_genomio_json_data_dir: Path, filename: str, error_pattern: str
) -> None:
    """Test that malformed RepeatMasker `.out` rows raise errors.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        filename: Name of RepeatMasker `.out` file containing invalid data.
        error_pattern: Regex expected in the raised error message.

    """
    out_path = convert_to_genomio_json_data_dir / filename
    with pytest.raises(ValueError, match=error_pattern):
        repeatmasker.parse_output(out_path, None)


@pytest.mark.parametrize(
    ("filename", "warning_pattern"),
    [
        pytest.param(
            "repeatmasker_out/error_non_positive_repeat_end.out",
            r"Invalid repeat coordinates",
            id="Skips non-positive repeat end with warning",
        ),
        pytest.param(
            "repeatmasker_out/error_repeat_end_before_start.out",
            r"repeat_end < repeat_start",
            id="Skips repeat end before start with warning",
        ),
    ],
)
def test_parse_output_skips_invalid_repeat_coordinates(
    convert_to_genomio_json_data_dir: Path,
    filename: str,
    warning_pattern: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test ``repeatmasker.parse_output()`` logs a warning and skips invalid records.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        filename: Input RepeatMasker `.out` filename containing invalid coordinates.
        warning_pattern: Substring expected to appear in the logged warning message.
        caplog: Pytest fixture for capturing log output.

    """
    output_path = convert_to_genomio_json_data_dir / filename

    with caplog.at_level("WARNING"):
        features, consensuses_by_key = repeatmasker.parse_output(output_path, None)

    assert not features
    assert not consensuses_by_key
    assert caplog.text
    assert any(warning_pattern in record.message for record in caplog.records)


@pytest.mark.parametrize(
    ("filename", "expected_fragments"),
    [
        pytest.param(
            "repeatmasker_out/error_multiple_errors.out",
            [
                "Found 3 errors while parsing RepeatMasker output",
                "Invalid 'score'",
                "seq_region_end < seq_region_start",
                "Unexpected strand token",
            ],
            id="RepeatMasker multiple row errors",
        ),
    ],
)
def test_parse_output_collates_all_errors(
    convert_to_genomio_json_data_dir: Path,
    filename: str,
    expected_fragments: list[str],
) -> None:
    """Test that malformed RepeatMasker parser rows are collated before raising.

    Args:
        convert_to_genomio_json_data_dir: Test data directory for convert-to-GenomIO JSON fixtures.
        filename: Input filename containing invalid data.
        expected_fragments: Expected substrings in the raised error message.

    """
    input_path = convert_to_genomio_json_data_dir / filename

    with pytest.raises(
        ValueError,
        match=r"^Found \d+ errors while parsing RepeatMasker output in .*:",
    ) as excinfo:
        repeatmasker.parse_output(input_path, None)

    error_message = str(excinfo.value)
    for expected_fragment in expected_fragments:
        assert expected_fragment in error_message
