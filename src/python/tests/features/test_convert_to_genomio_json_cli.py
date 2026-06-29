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
"""Unit testing of the convert-to-GenomIO JSON command-line entry point."""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json


@patch("ensembl.io.genomio.features.convert_to_genomio_json.cli.create_genomio_json")
@pytest.mark.parametrize(
    "use_repeatmodeler_lib",
    [
        pytest.param(False, id="Required arguments only, no consensus library"),
        pytest.param(True, id="All optional arguments overridden, consensus library used"),
    ],
)
def test_main_passes_expected_arguments(
    mock_create_genomio_json: Mock,
    convert_to_genomio_json_data_dir: Path,
    tmp_path: Path,
    use_repeatmodeler_lib: bool,
) -> None:
    """Test ``convert_to_genomio_json.main()`` passes parsed arguments correctly.

    Args:
        mock_create_genomio_json: Mock for ``convert_to_genomio_json.create_genomio_json()``.
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

    expected = {
        "input_path": input_path,
        "repeatmasker_consensus_lib_path": None,
        "output_path": output_path,
        "analysis_logic_name": "repeatmask_customlib",
        "analysis_display_label": "Repeats: Custom library",
        "analysis_description": (
            'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker'
            "</a>, using a custom library of <em>ab initio</em> repeat profiles for this species."
        ),
        "program": "RepeatMasker",
        "program_version": "4.1.5",
        "program_parameters": None,
        "source_provider": "Ensembl",
        "is_primary": False,
    }

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
        expected = {
            "input_path": input_path,
            "repeatmasker_consensus_lib_path": consensus_lib,
            "output_path": output_path,
            "analysis_logic_name": "repeatmask_repbase",
            "analysis_display_label": "Repeats: Repbase",
            "analysis_description": (
                'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker'
                '</a>, using the <a rel="external" href="http://www.girinst.org/repbase/">Repbase</a> '
                "library of repeat profiles."
            ),
            "program": "RepeatMasker",
            "program_version": "4.1.7",
            "program_parameters": "-lib foo",
            "source_provider": "Custom",
            "is_primary": True,
        }

    convert_to_genomio_json.main(argv)

    mock_create_genomio_json.assert_called_once_with(**expected)


@patch("ensembl.io.genomio.features.convert_to_genomio_json.cli.create_genomio_json")
def test_main_reraises_exceptions(
    mock_create_genomio_json: Mock, convert_to_genomio_json_data_dir: Path, tmp_path: Path
) -> None:
    """Test the ``convert_to_genomio_json.main()`` function reraises exceptions.

    Args:
        mock_create_genomio_json: Mock for ``convert_to_genomio_json.create_genomio_json()``.
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    output_path = tmp_path / "out.json"

    mock_create_genomio_json.side_effect = RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        convert_to_genomio_json.main(
            [
                "repeatmasker",
                "custom",
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--program-version",
                "4.1.5",
            ]
        )
