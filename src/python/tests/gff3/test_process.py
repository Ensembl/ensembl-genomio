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
"""Unit testing of `ensembl.io.genomio.gff3.process` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_process.py

"""

import filecmp
from pathlib import Path

import pytest

from ensembl.io.genomio.gff3.process import GFFGeneMerger


@pytest.mark.parametrize(
    "input_file, expected_file",
    [
        ("input1.gff3", "output1.gff3"),
    ],
)
def test_merge(data_dir: Path, tmp_path: Path, input_file: str, expected_file: str) -> None:
    """Tests the `GFFGeneMerger.merge()` method.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Fixture which will provide a temporary directory unique to the test invocation.
        input_file: Name of the GFF file with example input.
        expected_file: Name of the GFF file with the expected output.

    """
    merger = GFFGeneMerger()
    gff_input_path = data_dir / input_file
    test_output_path = tmp_path / f"{input_file}.test.gff3"
    expected_path = data_dir / expected_file
    merger.merge(gff_input_path, test_output_path)
    assert filecmp.cmp(test_output_path, expected_path)
