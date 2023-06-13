#!/usr/bin/env python
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
"""Unit testing of :mod:`ensembl.io.genomio.gff3.process_gff3` module.

"""

from pathlib import Path

import pytest

from ensembl.io.genomio.gff3.process_gff3 import GFFGeneMerger


class TestGFF3Merge:
    """Tests for the GFF3GeneMerger module."""

    tmp_dir: Path
    data_dir: Path

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, tmp_dir: Path, files_dir: Path):
        """Loads necessary fixtures and values as class attributes."""
        type(self).tmp_dir = tmp_dir
        type(self).data_dir = files_dir / "process_gff3/merge_genes"

    @pytest.mark.parametrize(
        "input_name, expected_name",
        [
            ("input1.gff3", "output1.gff3"),
        ],
    )
    def test_merge(self, tmp_path: Path, input_name: str, expected_name: str) -> None:
        """Tests merge method.

        Args:
            tmp_path: Where temporary files will be created.
            input_name: Name of the GFF file with example input, in the test folder.
            expected_name: Name of the GFF file with expected output, in the test folder.

        """
        merger = GFFGeneMerger()
        assert isinstance(merger, GFFGeneMerger)

        gff_input_path = self.data_dir / input_name
        test_output_path = tmp_path / f"{input_name}.test.gff3"

        merger.merge(gff_input_path, test_output_path)

        expected_path = self.data_dir / expected_name
        assert test_output_path.read_text() == expected_path.read_text()
