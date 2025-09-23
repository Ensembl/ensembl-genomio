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
"""Unit testing of `ensembl.io.genomio.seq_region.prepare` module."""
# pylint: disable=too-many-positional-arguments

from pathlib import Path
from typing import Callable
from unittest.mock import Mock, patch

import pytest

from ensembl.io.genomio.seq_region.prepare import prepare_seq_region_metadata


@pytest.mark.parametrize(
    ("gbff_path", "expected_path", "to_exclude"),
    [
        pytest.param(None, "out/no_gbff.json", [], id="Prepare without gbff"),
        pytest.param("apicoplast.gb", "out/with_gbff.json", [], id="Prepare with gbff"),
        pytest.param(None, "out/removed.json", ["NC_001799.1"], id="Prepare and exclude"),
    ],
)
@patch("ensembl.io.genomio.seq_region.collection.requests.get")
def test_prepare_seq_region_metadata(
    mock_requests_get: Mock,
    mock_response: Callable,
    data_dir: Path,
    tmp_path: Path,
    assert_files: Callable,
    gbff_path: Path | None,
    expected_path: Path,
    to_exclude: list[str],
) -> None:
    """Test `prepare_seq_region_metadata`.

    Args:
        mock_requests_get: Mock for requests.get.
        mock_response: Mock response for requests.
        data_dir: Module's test data directory fixture.
        tmp_path: Fixture which will provide a temporary directory unique to the test invocation.
        assert_files: Function to compare if two files are equal or not.
        gbff_path: Input GBFF file in any.
        expected_path: Expect JSON output.
        to_exclude: List of sequences to exclude.

    """
    gbff_file = None
    if gbff_path:
        gbff_file = data_dir / gbff_path
    dst_file = tmp_path / "output.json"
    mock_requests_get.return_value = mock_response("{}")  # Mock requests just in case
    prepare_seq_region_metadata(
        genome_file=data_dir / "genome.json",
        report_file=data_dir / "report.txt",
        gbff_file=gbff_file,
        mock_run=True,
        to_exclude=to_exclude,
        dst_file=dst_file,
    )
    assert_files(dst_file, data_dir / expected_path)
