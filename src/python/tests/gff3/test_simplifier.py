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
"""Unit testing of `ensembl.io.genomio.gff3.restructure` module."""

from contextlib import nullcontext as does_not_raise
from os import PathLike
from pathlib import Path
from typing import Any, Callable, ContextManager, Dict, List, Union

# from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
from pytest import param, raises

# from ensembl.io.genomio.gff3.exceptions import GFFParserError
from ensembl.io.genomio.gff3.simplifier import GFFSimplifier


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param("ok_genes.gff", "ok_genes_simped.gff", does_not_raise(), id="ok gene"),
    ],
)
def test_simpler_gff3(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    expectation: ContextManager,
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    with expectation:
        simp = GFFSimplifier()
        simp.simpler_gff3(input_gff)
        simp.records.to_gff(output_gff)
        assert_files(output_gff, data_dir / expected_gff)
