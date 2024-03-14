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
from typing import Callable, ContextManager, Dict, Optional

from Bio.SeqFeature import SeqFeature
import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.exceptions import GFFParserError
from ensembl.io.genomio.gff3.simplifier import GFFSimplifier


def check_one_feature(input_gff: PathLike, output_gff: PathLike, check_function: str) -> SeqFeature:
    """Load 1 feature from a GFF, apply a function, then write it back to a GFF."""
    simp = GFFSimplifier()
    simp.records.from_gff(input_gff)
    # Get the only feature
    feat = simp.records[0].features[0]
    # Apply the named function
    check_method = getattr(simp, check_function)
    new_feat = check_method(feat)

    if new_feat is not None:
        # Put it back
        simp.records[0].features = [new_feat]
        simp.records.to_gff(output_gff)
    return new_feat


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param("ok_gene.gff", "ok_gene.gff", does_not_raise(), id="ok gene"),
        param("lone_transcript.gff", "lone_transcript_simped.gff", does_not_raise(), id="lone transcript"),
        param("lone_trna.gff", "lone_trna_simped.gff", does_not_raise(), id="lone tRNA"),
        param("lone_rrna.gff", "lone_rrna_simped.gff", does_not_raise(), id="lone rRNA"),
    ],
)
def test_create_gene_for_lone_transcript(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    expectation: ContextManager,
) -> None:
    """Test gene create gene for lone transcript."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    with expectation:
        new_feat = check_one_feature(input_gff, output_gff, "create_gene_for_lone_transcript")
        if new_feat:
            assert_files(output_gff, Path(data_dir / expected_gff))


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param("ok_gene.gff", "ok_gene.gff", does_not_raise(), id="ok gene"),
        param("lone_cds.gff", "lone_cds_simped.gff", does_not_raise(), id="lone CDS"),
        param("lone_cds_pseudo.gff", "lone_cds_pseudo_simped.gff", does_not_raise(), id="lone pseudo CDS"),
    ],
)
def test_create_gene_for_lone_cds(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    expectation: ContextManager,
) -> None:
    """Test gene create gene for lone CDS."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    with expectation:
        new_feat = check_one_feature(input_gff, output_gff, "create_gene_for_lone_cds")
        if new_feat:
            assert_files(output_gff, Path(data_dir / expected_gff))


@pytest.mark.parametrize(
    "in_type, in_mobile_type, in_product, out_type, out_description, expectation",
    [
        param("gene", None, None, "gene", None, does_not_raise(), id="Gene, skip"),
        param("transposable_element", None, None, "transposable_element", None, does_not_raise(), id="TE"),
        param("mobile_genetic_element", None, None, "transposable_element", None, does_not_raise(), id="MGE"),
        param("transposable_element", "transposon", None, "transposable_element", "transposon", does_not_raise(), id="MGE, transposon"),
        param("transposable_element", "transposon:LOREM", None, "transposable_element", "transposon (LOREM)", does_not_raise(), id="MGE, transposon named"),
        param("transposable_element", "retrotransposon:LOREM", None, "transposable_element", "retrotransposon (LOREM)", does_not_raise(), id="MGE, retrotransposon named"),
        param("transposable_element", "UNKNOWNtransposon:LOREM", None, "transposable_element", None, raises(GFFParserError), id="MGE, unknown type"),
        param("transposable_element", "transposon", "PROD", "transposable_element", "PROD", does_not_raise(), id="MGE, transposon, product exists"),
    ],
)
def test_normalize_non_gene(
    in_type: str,
    in_mobile_type: Optional[str],
    in_product: Optional[str],
    out_type: str,
    out_description: Optional[str],
    expectation: ContextManager,
) -> None:
    """Test gene create gene for lone CDS."""
    simp = GFFSimplifier()
    feat = SeqFeature(None, in_type)
    feat.qualifiers = {"source": "LOREM"}
    if in_mobile_type is not None:
        feat.qualifiers["mobile_element_type"] = [in_mobile_type]
    if in_product is not None:
        feat.qualifiers["product"] = [in_product]
    feat.sub_features = []
    with expectation:
        new_feat = simp.normalize_non_gene(feat)
        if new_feat is not None:
            assert new_feat.type == out_type
            if out_description is not None:
                assert new_feat.qualifiers == feat.qualifiers


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param("ok_gene.gff", "ok_gene.gff", does_not_raise(), id="ok gene"),
        param("gene_ignored.gff", None, does_not_raise(), id="gene ignored"),
        param("mobile_te.gff", "mobile_te.gff", does_not_raise(), id="TE"),
        param("ok_protein_coding_gene.gff", "ok_gene.gff", does_not_raise(), id="ok protein_coding_gene"),
    ],
)
def test_simpler_gff3_feature(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: Optional[PathLike],
    expectation: ContextManager,
) -> None:
    """Test simplifying one gene (from a GFF3 file)."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    with expectation:
        new_feat = check_one_feature(input_gff, output_gff, "simpler_gff3_feature")
        if expected_gff is None:
            assert new_feat is None
        else:
            assert_files(output_gff, Path(data_dir / expected_gff))


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param("ok_gene.gff", "ok_gene.gff", does_not_raise(), id="ok gene"),
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
        param("ok_gene.gff", "ok_gene.gff", False, does_not_raise(), id="ok type, Keep"),
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
        param(
            None,
            "genes_badnames.gff",
            "genes_badnames_noname.gff",
            does_not_raise(),
            id="Genes with bad names, no genome",
        ),
        param(
            "genome_no_brc4.json",
            "genes_badnames.gff",
            "genes_badnames_noname.gff",
            does_not_raise(),
            id="Genes with bad names, genome not BRC4",
        ),
        param(
            "genome_brc4.json",
            "genes_badnames.gff",
            "genes_badnames_brc4name.gff",
            does_not_raise(),
            id="Genes with bad names, genome BRC4",
        ),
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
