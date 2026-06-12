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
"""Unit testing of `ensembl.io.genomio.metadata` module."""

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager
from unittest.mock import MagicMock, patch

from deepdiff import DeepDiff
import pytest
import yaml

from ensembl.io.genomio import metadata


SAMPLE_ACRONYMS = {
    "GenBank": "insdc",
    "RefSeq": "insdc",
    "Ensembl": "ebi",
}

EXPECTED_YAML_CONTENT = {
    "insdc": ["GenBank", "RefSeq"],
    "ebi": ["Ensembl"],
}


@patch("ensembl.io.genomio.metadata.get_json", return_value=SAMPLE_ACRONYMS)
def test_create_provider_ftp_yaml(mock_get_json: MagicMock, tmp_path: Path) -> None:
    """Test :func:`metadata.create_provider_ftp_yaml`.

    Args:
        mock_get_json: Patched `get_json()` returning a controlled mapping.
        tmp_path: Pytest-provided temporary directory.

    """
    output = tmp_path / "providers.yaml"
    metadata.create_provider_ftp_yaml(output)
    assert output.exists(), "Output YAML file was not created"
    with output.open("r") as fh:
        content = yaml.safe_load(fh)
    assert not DeepDiff(content, EXPECTED_YAML_CONTENT), "YAML content does not match expected structure"


@patch("ensembl.io.genomio.metadata.get_json", return_value=SAMPLE_ACRONYMS)
def test_create_provider_ftp_yaml_accepts_str(mock_get_json: MagicMock, tmp_path: Path) -> None:
    """Test :func:`metadata.create_provider_ftp_yaml` accepts plain string paths.

    Args:
        mock_get_json: Patched `get_json()` returning a controlled mapping.
        tmp_path: Pytest-provided temporary directory.

    """
    output = tmp_path / "providers.yaml"
    metadata.create_provider_ftp_yaml(str(output))
    assert output.exists(), "Output YAML file was not created"


@pytest.mark.parametrize(
    ("argv", "expectation"),
    [
        pytest.param(
            ["--output", "out.yaml"],
            does_not_raise({"output": "out.yaml", "log_level": "WARNING"}),
            id="Default args"
        ),
        pytest.param([], pytest.raises(SystemExit), id="Missing required args"),
    ],
)
def test_parse_args(argv: list[str] | None, expectation: ContextManager) -> None:
    """Test :func:`metadata.parse_args`.

    Args:
        argv: Argument vector to parse.
        expectation: Context manager for the expected outcome.

    """
    with expectation as expected:
        args = metadata.parse_args(argv)
        # DeepDiff is not able to compare two objects of Path type, so convert it to string
        args.output = str(args.output)
        assert not DeepDiff(vars(args), expected)


@patch("ensembl.io.genomio.metadata.create_provider_ftp_yaml")
def test_main(mock_create: MagicMock) -> None:
    """Test :func:`metadata.main`.

    Args:
        mock_create: Patched `create_provider_ftp_yaml()` to verify it is called with the expected
            output path.

    """
    metadata.main(["--output", "out.yaml"])
    mock_create.assert_called_once_with(output_path=Path("out.yaml"))


@patch("ensembl.io.genomio.metadata.create_provider_ftp_yaml")
def test_main_raises_exception(mock_create: MagicMock) -> None:
    """Test :func:`metadata.main`.

    Args:
        mock_create: Patched `create_provider_ftp_yaml()` to verify it is called with the expected
            output path.

    """
    mock_create.side_effect = RuntimeError("Mocked exception")
    with pytest.raises(RuntimeError, match=r"Mocked exception"):
        metadata.main(["--output", "out.yaml"])
