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

"""Tests for the repeatmasker_out_to_genomio_json module."""

import hashlib
import json
from pathlib import Path
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
    norm_name = name.strip()
    norm_class = repeat_class.strip()
    norm_type = repeat_type.strip()
    norm_seq = "".join(seq.split()).upper()
    payload = f"{norm_name}\t{norm_class}\t{norm_type}\t{norm_seq}".encode("utf-8")
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
    ("repeat_type", "expected"),
    [
        ("LINE1", "Type I Transposons/LINE"),
        ("Satellite", "Satellite repeats"),
        ("fooRNA", "RNA repeats"),
        ("NotInMap", "Unknown"),
    ],
)
def test_map_rm_repeat_type(repeat_type: str, expected: str) -> None:
    """
    Tests the ``convert_to_genomio_json._map_rm_repeat_type()`` function for known and unknown patterns.

    Args:
        repeat_type: Input repeat type string.
        expected: Expected mapped output.
    """
    assert convert_to_genomio_json._map_rm_repeat_type(repeat_type) == expected


def test_parse_rm_consensus_library(data_dir: Path) -> None:
    """
    Tests the ``convert_to_genomio_json._parse_rm_consensus_library()`` internal helper.

    Args:
        data_dir: Module's test data directory fixture.
    """
    fasta_path = data_dir / "repeatmodeler" / "classified_mixed.fa"

    expected_consensuses = {
        _sha256_key("AluY", "SINE", "Unknown", "acgt"): convert_to_genomio_json.Consensus(
            name="AluY",
            repeat_class="SINE",
            repeat_type="Unknown",
            seq="acgt",
        ),
        _sha256_key("Foo", "LINE", "Unknown", "AATT"): convert_to_genomio_json.Consensus(
            name="Foo",
            repeat_class="LINE",
            repeat_type="Unknown",
            seq="AATT",
        ),
        _sha256_key("Bar", "Unknown", "Unknown", "ccgg"): convert_to_genomio_json.Consensus(
            name="Bar",
            repeat_class="Unknown",
            repeat_type="Unknown",
            seq="ccgg",
        ),
    }

    expected_keys_by_triplet = {
        ("AluY", "SINE", "Alu"): _sha256_key("AluY", "SINE", "Unknown", "acgt"),
        ("Foo", "LINE", "Unknown"): _sha256_key("Foo", "LINE", "Unknown", "AATT"),
        ("Bar", "Unknown", "Unknown"): _sha256_key("Bar", "Unknown", "Unknown", "ccgg"),
    }

    consensus_keys_by_triplet, consensuses_by_key = convert_to_genomio_json._parse_rm_consensus_library(
        fasta_path
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
                        "rm_repeat_type": "Alu",
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
                        "rm_repeat_type": "L1",
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
                        "rm_repeat_type": "Unknown",
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
                        "rm_repeat_type": "Alu",
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
                        "rm_repeat_type": "Alu",
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
                        "rm_repeat_type": "Unknown",
                    },
                },
            ],
            id="Skips blank lines",
        ),
    ],
)
def test_parse_repeatmasker_out_success(
    data_dir: Path,
    filename: str,
    expected_features: list[dict[str, object]],
) -> None:
    """
    Tests successful parsing of RepeatMasker ``.out`` files.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Input RepeatMasker ``.out`` filename.
        expected_features: Expected parsed repeat feature dictionaries.
    """
    out_path = data_dir / filename
    features, consensuses_by_key = convert_to_genomio_json.parse_repeatmasker_out(
        out_path,
        rm_consensus_lib_path=None,
    )

    assert features == expected_features
    assert consensuses_by_key == {}


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
def test_parse_repeatmasker_out_errors(data_dir: Path, filename: str, error_pattern: str) -> None:
    """
    Tests that malformed RepeatMasker ``.out`` rows raise errors.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Name of RepeatMasker ``.out`` file containing invalid data.
        error_pattern: Regex expected in the raised error message.
    """
    out_path = data_dir / filename
    with pytest.raises(ValueError, match=error_pattern):
        convert_to_genomio_json.parse_repeatmasker_out(out_path, None)


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
def test_parse_repeatmasker_out_skips_invalid_repeat_coordinates(
    data_dir: Path,
    filename: str,
    warning_pattern: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """
    Tests that `convert_to_genomio_json.parse_repeatmasker_out()` logs a warning and skips invalid records.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Input RepeatMasker ``.out`` filename containing invalid coordinates.
        warning_pattern: Substring expected to appear in the logged warning message.
        caplog: Pytest fixture for capturing log output.
    """
    out_path = data_dir / filename

    with caplog.at_level("WARNING"):
        features, consensuses_by_key = convert_to_genomio_json.parse_repeatmasker_out(out_path, None)

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
def test_parse_trf_out_success(
    data_dir: Path,
    filename: str,
    expected_features: list[dict[str, object]],
    expected_consensuses: dict[str, convert_to_genomio_json.Consensus],
) -> None:
    """
    Tests successful parsing of TRF ``.dat`` examples.

    Args:
        data_dir: Module's test data directory fixture.
        filename: Input TRF ``.dat`` filename.
        expected_features: Expected parsed repeat feature dictionaries.
        expected_consensuses: Expected parsed consensus records.
    """
    trf_path = data_dir / filename

    features, consensuses_by_key = convert_to_genomio_json.parse_trf_out(trf_path)

    assert features == expected_features
    assert consensuses_by_key == expected_consensuses


def test_parse_trf_out_errors_without_sequence_header(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Tests ``convert_to_genomio_json.parse_trf_out()`` raises when lines are parsed without a sequence header.

    Args:
        tmp_path: Temporary directory provided by pytest.
        monkeypatch: Pytest fixture for temporarily patching module attributes.
    """

    class FakeSeqMatch:
        @staticmethod
        def group(name: str) -> str | None:
            return None

    class FakeSequenceRegex:
        @staticmethod
        def match(line: str) -> FakeSeqMatch | None:
            if line == "Sequence:":
                return FakeSeqMatch()
            return None

    trf_path = tmp_path / "missing_sequence_name.dat"
    trf_path.write_text(
        "Sequence:\n1 4 4 1.0 4 100 0 50 25 25 25 25 2.0 AT\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(convert_to_genomio_json, "TRF_SEQUENCE_RE", FakeSequenceRegex)

    with pytest.raises(ValueError, match=r"TRF sequence header not found before data lines"):
        convert_to_genomio_json.parse_trf_out(trf_path)


@pytest.mark.parametrize(
    (
        "input_filename",
        "analysis_logic_name",
        "rm_consensus_lib_filename",
        "expected_repeat_consensus",
        "expected_feature",
    ),
    [
        param(
            "create_json/basic.out",
            "repeatmask_customlib",
            None,
            None,
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
                    "rm_repeat_type": "Alu",
                },
            },
            id="RepeatMasker without RepeatModeler consensuslibrary",
        ),
        param(
            "create_json/basic.out",
            "repeatmask_customlib",
            "create_json/repeatmodeler_match.fa",
            [
                {
                    "repeat_consensus_key": _sha256_key("Foo", "SINE", "Unknown", "ACGT"),
                    "repeat_name": "Foo",
                    "repeat_class": "SINE",
                    "repeat_type": "Unknown",
                    "repeat_consensus": "acgt",
                }
            ],
            {
                "seq_region": "chr1",
                "seq_region_start": 10,
                "seq_region_end": 20,
                "seq_region_strand": "+",
                "repeat_start": 1,
                "repeat_end": 11,
                "repeat_consensus": _sha256_key("Foo", "SINE", "Unknown", "ACGT"),
                "score": 100.0,
                "attributes": {
                    "perc_div": 1.0,
                    "perc_del": 0.0,
                    "perc_ins": 0.0,
                    "rm_repeat_type": "Alu",
                },
            },
            id="RepeatMasker with RepeatModeler consensus library",
        ),
        param(
            "trf/success_window_with_params.dat",
            "trf",
            None,
            [
                {
                    "repeat_consensus_key": _sha256_key("trf", "trf", "Tandem repeats", "AT"),
                    "repeat_name": "trf",
                    "repeat_class": "trf",
                    "repeat_type": "Tandem repeats",
                    "repeat_consensus": "AT",
                }
            ],
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
            id="Tandem Repeat Finder",
        ),
        param(
            "create_json/basic.out",
            "unsupported",
            None,
            None,
            None,
            id="Unsupported logic name",
        ),
    ],
)
def test_create_genomio_json(
    data_dir: Path,
    tmp_path: Path,
    input_filename: str,
    analysis_logic_name: str,
    rm_consensus_lib_filename: str | None,
    expected_repeat_consensus: list[dict[str, str]] | None,
    expected_feature: dict[str, object],
) -> None:
    """
    Tests JSON creation with and without consensus library input.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        input_filename: Input feature filename.
        analysis_logic_name: Analysis logic name.
        rm_consensus_lib_filename: Optional RepeatModeler consensus library filename.
        expected_repeat_consensus: Expected repeat consensus JSON section.
        expected_feature: Expected feature JSON section.
    """
    input = data_dir / input_filename
    rm_consensus_lib = data_dir / rm_consensus_lib_filename if rm_consensus_lib_filename is not None else None
    output = tmp_path / "out.json"

    convert_to_genomio_json.create_genomio_json(
        input_path=input,
        rm_consensus_lib_path=rm_consensus_lib,
        output_path=output,
        analysis_logic_name=analysis_logic_name,
        analysis_display_label="Repeats: Custom library",
        analysis_description="desc",
        program="RepeatMasker",
        program_version="4.1.5",
        program_parameters="-nolow -gccalc -q",
        source_provider="Ensembl",
        is_primary=False,
    )

    doc = json.loads(output.read_text(encoding="utf-8"))

    assert doc["analysis"]["logic_name"] == analysis_logic_name
    assert doc["analysis"]["display_label"] == "Repeats: Custom library"
    assert doc["analysis"]["description"] == "desc"
    assert doc["analysis"]["program"] == "RepeatMasker"
    assert doc["analysis"]["program_version"] == "4.1.5"
    assert doc["analysis"]["program_parameters"] == "-nolow -gccalc -q"
    assert doc["analysis"]["run_date"].endswith("Z")

    assert doc["source"] == {"source_provider": "Ensembl", "is_primary": False}

    if expected_feature is None:
        assert doc["repeat_features"] == []
    else:
        assert doc["repeat_features"] == [expected_feature]

    if expected_repeat_consensus is None:
        assert "repeat_consensus" not in doc
    else:
        assert doc["repeat_consensus"] == expected_repeat_consensus


def test_create_genomio_json_omits_program_parameters_when_none(data_dir: Path, tmp_path: Path) -> None:
    """
    Tests that program parameters are omitted when no value is provided.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
    """
    input = data_dir / "create_json" / "basic_dna.out"
    output = tmp_path / "out.json"

    convert_to_genomio_json.create_genomio_json(
        input_path=input,
        rm_consensus_lib_path=None,
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
    Tests the ``convert_to_genomio_json.parse_args()`` function with real input files.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        use_repeatmodeler_lib: Whether to include optional repeatmodeler/program-parameters args.
    """
    input = data_dir / "create_json" / "basic.out"
    output = tmp_path / "out.json"

    argv = [
        "--input",
        str(input),
        "--output",
        str(output),
        "--analysis-logic-name",
        "repeatmask_customlib",
        "--program-version",
        "4.1.5",
    ]

    if use_repeatmodeler_lib:
        rm_consensus_lib = data_dir / "create_json" / "repeatmodeler_match.fa"
        argv = [
            "--input",
            str(input),
            "--rm-consensus-lib",
            str(rm_consensus_lib),
            "--output",
            str(output),
            "--analysis-logic-name",
            "repeatmask_rmlib",
            "--analysis-display-label",
            "Repeats: RepBase",
            "--analysis-description",
            "RepBase description",
            "--program",
            "RepeatMasker",
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
        assert args.analysis_display_label == "Repeats: RepBase"
        assert args.analysis_description == "RepBase description"
        assert args.program == "RepeatMasker"
        assert args.source_provider == "Custom"
        assert args.is_primary is True
        assert args.rm_consensus_lib == data_dir / "create_json" / "repeatmodeler_match.fa"
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
        assert not hasattr(args, "rm_consensus_lib")
        assert not hasattr(args, "program_parameters")


@pytest.mark.parametrize(
    ("analysis_logic_name", "include_rm_consensus_lib", "error_pattern"),
    [
        param(
            "unsupported",
            False,
            r"Unsupported analysis logic name",
            id="Unsupported analysis logic name",
        ),
        param(
            "trf",
            True,
            r"--rm-consensus-lib is not applicable",
            id="Rejects --rm-consensus-lib with non-RepeatMasker logic name",
        ),
    ],
)
def test_parse_args_errors(
    data_dir: Path,
    tmp_path: Path,
    analysis_logic_name: str,
    include_rm_consensus_lib: bool,
    error_pattern: str,
) -> None:
    """
    Tests ``convert_to_genomio_json.parse_args()`` rejects unsupported argument combinations.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        analysis_logic_name: Analysis logic name to parse.
        include_rm_consensus_lib: Whether to include a RepeatMasker consensus library path.
        error_pattern: Regex expected in the raised error message.
    """
    input = data_dir / "create_json" / "basic.out"
    output = tmp_path / "out.json"
    argv = [
        "--input",
        str(input),
        "--output",
        str(output),
        "--analysis-logic-name",
        analysis_logic_name,
        "--program-version",
        "4.1.5",
    ]
    if include_rm_consensus_lib:
        argv.extend(
            [
                "--rm-consensus-lib",
                str(data_dir / "create_json" / "repeatmodeler_match.fa"),
            ]
        )

    with pytest.raises(ValueError, match=error_pattern):
        convert_to_genomio_json.parse_args(argv)


@patch("ensembl.io.genomio.features.convert_to_genomio_json.create_genomio_json")
@pytest.mark.parametrize(
    "use_repeatmodeler_lib",
    [
        param(False, id="Required arguments only"),
        param(True, id="All optional arguments overridden"),
    ],
)
def test_main_passes_expected_arguments(
    mock_create_genomio_json: Mock,
    data_dir: Path,
    tmp_path: Path,
    use_repeatmodeler_lib: bool,
) -> None:
    """
    Tests ``convert_to_genomio_json.main()`` passes parsed arguments through to `create_genomio_json()` using real input files.

    Args:
        mock_create_genomio_json: Mock for `create_genomio_json()`.
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        use_repeatmodeler_lib: Whether to include optional repeatmodeler/program-parameters args.
    """
    input = data_dir / "create_json" / "basic.out"
    output = tmp_path / "out.json"

    argv = [
        "--input",
        str(input),
        "--output",
        str(output),
        "--analysis-logic-name",
        "repeatmask_customlib",
        "--program-version",
        "4.1.5",
    ]

    expected = {
        "input_path": input,
        "rm_consensus_lib_path": None,
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
        rm_consensus_lib = data_dir / "create_json" / "repeatmodeler_match.fa"
        argv = [
            "--input",
            str(input),
            "--rm-consensus-lib",
            str(rm_consensus_lib),
            "--output",
            str(output),
            "--analysis-logic-name",
            "repeatmask_rmlib",
            "--analysis-display-label",
            "RepBase",
            "--analysis-description",
            "RepBase desc",
            "--program",
            "RepeatMasker",
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
            "rm_consensus_lib_path": rm_consensus_lib,
            "output_path": output,
            "analysis_logic_name": "repeatmask_rmlib",
            "analysis_display_label": "RepBase",
            "analysis_description": "RepBase desc",
            "program": "RepeatMasker",
            "program_version": "4.1.7",
            "program_parameters": "-lib foo",
            "source_provider": "Custom",
            "is_primary": True,
        }

    assert convert_to_genomio_json.main(argv) is None
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
                "--input",
                str(input),
                "--output",
                str(output),
                "--analysis-logic-name",
                "repeatmask_customlib",
                "--program-version",
                "4.1.5",
            ]
        )
