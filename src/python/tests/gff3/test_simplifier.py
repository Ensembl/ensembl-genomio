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
from typing import Any, Callable, ContextManager, Dict, List, Optional, Union

# from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.exceptions import GFFParserError
from ensembl.io.genomio.gff3.simplifier import GFFSimplifier


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param("ok_genes.gff", "ok_genes_simped.gff", does_not_raise(), id="ok gene"),
        param("bad_gene_type.gff", "", raises(GFFParserError), id="Unsupported gene type"),
        param("bad_tr_type.gff", "", raises(GFFParserError), id="Unsupported transcript type"),
        param("bad_subtr_type.gff", "", raises(GFFParserError), id="Unsupported subtranscript type"),
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


@pytest.mark.parametrize(
    "in_gff, expected_gff, skip_unrecognized, expectation",
    [
        param("bad_gene_type.gff", "", None, raises(GFFParserError), id="Unset skip unrecognized, fail"),
        param("bad_gene_type.gff", "", True, raises(GFFParserError), id="True skip unrecognized, fail"),
        param("bad_gene_type.gff", "bad_gene_type.gff", False, does_not_raise(), id="bad type, Keep"),
        param("ok_genes.gff", "ok_genes_simped.gff", False, does_not_raise(), id="ok type, Keep"),
    ],
)
def test_simpler_gff3_skip(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    skip_unrecognized: Optional[bool],
    expectation: ContextManager,
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    with expectation:
        simp = GFFSimplifier()
        if skip_unrecognized is not None:
            simp.skip_unrecognized = skip_unrecognized
        simp.simpler_gff3(input_gff)
        simp.records.to_gff(output_gff)
        assert_files(output_gff, data_dir / expected_gff)

@pytest.mark.parametrize(
    "genome_file, in_gff, expected_gff, expectation",
    [
        param(None, "genes_badnames.gff", "genes_badnames_noname.gff", does_not_raise(), id="Genes with bad names, no genome"),
        param("genome_no_brc4.json", "genes_badnames.gff", "genes_badnames_noname.gff", does_not_raise(), id="Genes with bad names, genome not BRC4"),
        param("genome_brc4.json", "genes_badnames.gff", "genes_badnames_brc4name.gff", does_not_raise(), id="Genes with bad names, genome BRC4"),
    ],
)
def test_gffsimplifier_with_genome(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    genome_file: Optional[PathLike],
    in_gff: PathLike,
    expected_gff: PathLike,
    expectation: ContextManager,
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    with expectation:
        if genome_file is None:
            simp = GFFSimplifier()
        else:
            simp = GFFSimplifier(genome_path=data_dir / genome_file)
        simp.simpler_gff3(input_gff)
        simp.records.to_gff(output_gff)
        assert_files(output_gff, data_dir / expected_gff)
