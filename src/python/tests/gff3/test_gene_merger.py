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
"""Unit testing of `ensembl.io.genomio.gff3.gene_merger` module."""

from pathlib import Path
from typing import Callable

import pytest

from ensembl.io.genomio.gff3 import GFFGeneMerger


@pytest.mark.parametrize(
    "input_file, expected_file",
    [
        pytest.param("merge_split_in.gff3", "merge_split_out.gff3", id="Split gene"),
        pytest.param("merge_fasta_in.gff3", "merge_fasta_out.gff3", id="Split gene with fasta"),
        pytest.param("merge_gene_only_in.gff3", "merge_gene_only_out.gff3", id="Split gene only"),
        pytest.param("merge_unordered_in.gff3", "merge_unordered_out.gff3", id="Unordered split gene"),
    ],
)
def test_merge(
    assert_files: Callable[[Path, Path], None],
    data_dir: Path,
    tmp_path: Path,
    input_file: str,
    expected_file: str,
) -> None:
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
    merger.merge(gff_input_path, test_output_path)

    expected_path = data_dir / expected_file
    assert_files(test_output_path, expected_path)
