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
"""Unit testing of convert-to-GenomIO JSON argument parsing."""

from pathlib import Path

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json


@pytest.mark.parametrize(
    "use_repeatmodeler_lib",
    [
        pytest.param(False, id="Required arguments only"),
        pytest.param(True, id="All optional arguments overridden"),
    ],
)
def test_parse_args(
    convert_to_genomio_json_data_dir: Path, tmp_path: Path, use_repeatmodeler_lib: bool
) -> None:
    """Test the ``convert_to_genomio_json.parse_args()`` function.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        use_repeatmodeler_lib: Whether to include optional repeatmodeler/program-parameters args.

    """
    input_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    output_path = tmp_path / "out.json"

    argv = [
        "repeatmasker",
        "custom",
        "--input",
        str(input_path),
        "--output",
        str(output_path),
        "--program-version",
        "4.1.5",
    ]

    if use_repeatmodeler_lib:
        consensus_lib = convert_to_genomio_json_data_dir / "create_json" / "repeatmodeler_match.fa"
        argv = [
            "repeatmasker",
            "repbase",
            "--input",
            str(input_path),
            "--consensus-lib",
            str(consensus_lib),
            "--output",
            str(output_path),
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
    assert args.input == input_path
    assert args.output == output_path

    if use_repeatmodeler_lib:
        assert args.program_version == "4.1.7"
        assert args.analysis_logic_name == "repeatmask_repbase"
        assert args.analysis_display_label == "Repeats: Repbase"
        assert args.analysis_description == (
            'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker'
            '</a>, using the <a rel="external" href="http://www.girinst.org/repbase/">Repbase</a> '
            "library of repeat profiles."
        )
        assert args.program == "RepeatMasker"
        assert args.source_provider == "Custom"
        assert args.is_primary is True
        assert (
            args.consensus_lib == convert_to_genomio_json_data_dir / "create_json" / "repeatmodeler_match.fa"
        )
        assert args.program_parameters == "-lib foo"
    else:
        assert args.program_version == "4.1.5"
        assert args.analysis_logic_name == "repeatmask_customlib"
        assert args.analysis_display_label == "Repeats: Custom library"
        assert args.analysis_description == (
            'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker'
            "</a>, using a custom library of <em>ab initio</em> repeat profiles for this species."
        )
        assert args.program == "RepeatMasker"
        assert args.source_provider == "Ensembl"
        assert args.is_primary is False
        assert not hasattr(args, "consensus_lib")
        assert not hasattr(args, "program_parameters")


def test_parse_args_red(convert_to_genomio_json_data_dir: Path, tmp_path: Path) -> None:
    """Test the ``convert_to_genomio_json.parse_args()`` function for Red input.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    output_path = tmp_path / "out.json"

    args = convert_to_genomio_json.parse_args(
        [
            "red",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--program-version",
            "2.0",
            "--program-parameters",
            "-gnm genome.fa",
            "--source-provider",
            "Custom",
            "--is-primary",
        ]
    )

    assert args.__class__.__name__ == "Namespace"
    assert args.input == input_path
    assert args.output == output_path
    assert args.program_version == "2.0"
    assert args.analysis_logic_name == "repeatdetector"
    assert args.analysis_display_label == "Repeats: Red"
    assert args.analysis_description == (
        'Repeats detected using <a href="https://bmcbioinformatics.biomedcentral.com/articles'
        '/10.1186/s12859-015-0654-5">Red (REPeatDetector)</a>'
    )
    assert args.program == "Red"
    assert args.repeatmasker_consensus_lib_path is None
    assert args.source_provider == "Custom"
    assert args.is_primary is True
    assert not hasattr(args, "consensus_lib")
    assert args.program_parameters == "-gnm genome.fa"


def test_parse_args_trf(convert_to_genomio_json_data_dir: Path, tmp_path: Path) -> None:
    """Test the ``convert_to_genomio_json.parse_args()`` function for TRF input.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    output_path = tmp_path / "out.json"

    args = convert_to_genomio_json.parse_args(
        [
            "trf",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--program-version",
            "4.10.0",
        ]
    )

    assert args.__class__.__name__ == "Namespace"
    assert args.input == input_path
    assert args.output == output_path
    assert args.program_version == "4.10.0"
    assert args.analysis_logic_name == "trf"
    assert args.analysis_display_label == "Tandem repeats (TRF)"
    assert args.analysis_description == (
        '<a rel="external" href="https://tandem.bu.edu/trf/trf.html">Tandem Repeats Finder</a> '
        "locates adjacent copies of a pattern of nucleotides."
    )
    assert args.program == "trf"
    assert args.repeatmasker_consensus_lib_path is None
    assert args.source_provider == "Ensembl"
    assert args.is_primary is False
    assert not hasattr(args, "consensus_lib")
    assert not hasattr(args, "program_parameters")
