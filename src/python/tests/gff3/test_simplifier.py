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
from typing import Callable, ContextManager, Optional

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
        if isinstance(new_feat, list):
            simp.records[0].features = new_feat
        else:
            simp.records[0].features = [new_feat]
        simp.records.to_gff(output_gff)
    return new_feat


@pytest.mark.parametrize(
    "in_gff, expected_gff",
    [
        param("ok_gene.gff", "ok_gene.gff", id="ok gene"),
        param("lone/transcript.gff", "lone/transcript_simped.gff", id="lone transcript"),
        param("lone/trna.gff", "lone/trna_simped.gff", id="lone tRNA"),
        param("lone/rrna.gff", "lone/rrna_simped.gff", id="lone rRNA"),
    ],
)
def test_create_gene_for_lone_transcript(
    tmp_path: Path,
    data_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
) -> None:
    """Test gene create gene for lone transcript."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / Path(in_gff).name
    new_feat = check_one_feature(input_gff, output_gff, "create_gene_for_lone_transcript")
    if new_feat:
        assert_files(output_gff, Path(data_dir / expected_gff))


@pytest.mark.parametrize(
    "in_gff, expected_gff",
    [
        param("ok_gene.gff", "ok_gene.gff", id="ok gene"),
        param("lone/cds.gff", "lone/cds_simped.gff", id="lone CDS"),
        param("lone/cds_pseudo.gff", "lone/cds_pseudo_simped.gff", id="lone pseudo CDS"),
    ],
)
def test_create_gene_for_lone_cds(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
) -> None:
    """Test gene create gene for lone CDS."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / Path(in_gff).name
    new_feat = check_one_feature(input_gff, output_gff, "create_gene_for_lone_cds")
    if new_feat:
        assert_files(output_gff, Path(data_dir / expected_gff))


@pytest.mark.parametrize(
    "in_type, in_mobile_type, in_product, out_type, out_description, expectation",
    [
        param("gene", None, None, "gene", None, does_not_raise(), id="Gene, skip"),
        param("transposable_element", None, None, "transposable_element", None, does_not_raise(), id="TE"),
        param("mobile_genetic_element", None, None, "transposable_element", None, does_not_raise(), id="MGE"),
        param(
            "transposable_element",
            "transposon",
            None,
            "transposable_element",
            "transposon",
            does_not_raise(),
            id="MGE, transposon",
        ),
        param(
            "transposable_element",
            "transposon:LOREM",
            None,
            "transposable_element",
            "transposon (LOREM)",
            does_not_raise(),
            id="MGE, transposon named",
        ),
        param(
            "transposable_element",
            "retrotransposon:LOREM",
            None,
            "transposable_element",
            "retrotransposon (LOREM)",
            does_not_raise(),
            id="MGE, retrotransposon named",
        ),
        param(
            "transposable_element",
            "UNKNOWNtransposon:LOREM",
            None,
            "transposable_element",
            None,
            raises(GFFParserError),
            id="MGE, unknown type",
        ),
        param(
            "transposable_element",
            "transposon",
            "PROD",
            "transposable_element",
            "PROD",
            does_not_raise(),
            id="MGE, transposon, product exists",
        ),
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
    """Test non-gene normalization."""
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


def test_normalize_non_gene_not_implemented() -> None:
    """Test non-gene not in the biotype list."""
    simp = GFFSimplifier()
    simp._biotypes = {"non_gene": {"supported": ["non_gene_name"]}}  # pylint: disable=protected-access
    feat = SeqFeature(None, "non_gene_name")
    with raises(NotImplementedError):
        simp.normalize_non_gene(feat)


@pytest.mark.parametrize(
    "in_type, name, out_type, expectation",
    [
        param("mRNA", "", "mRNA", does_not_raise(), id="mRNA no change"),
        param("C_gene_segment", "", "C_gene_segment", raises(GFFParserError), id="no standard name"),
        param("C_gene_segment", "immunoglobulin", "IG_C_gene", does_not_raise(), id="C immunoglobulin"),
        param("C_gene_segment", "ig", "IG_C_gene", does_not_raise(), id="C ig"),
        param("V_gene_segment", "t-cell", "TR_V_gene", does_not_raise(), id="V t-cell"),
        param("V_gene_segment", "T_cell", "TR_V_gene", does_not_raise(), id="V T_cell"),
        param("V_gene_segment", "Lorem Ipsum", "", raises(GFFParserError), id="V T_cell"),
    ],
)
def test_format_gene_segments(
    in_type: str,
    name: str,
    out_type: str,
    expectation: ContextManager,
) -> None:
    """Test gene create gene for lone CDS."""
    simp = GFFSimplifier()
    feat = SeqFeature(None, in_type)
    if name:
        feat.qualifiers["standard_name"] = [name]
    feat.sub_features = []
    with expectation:
        new_feat = simp.format_gene_segments(feat)
        assert new_feat.type == out_type


@pytest.mark.parametrize(
    "in_gff, expected_gff",
    [
        param("ok_gene.gff", "ok_gene.gff", id="ok gene"),
        param("clean/extra.gff", "clean/extra_clean.gff", id="ok gene with extra attribs"),
    ],
)
def test_clean_gene(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
) -> None:
    """Test clean gene."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / Path(in_gff).name
    new_feat = check_one_feature(input_gff, output_gff, "clean_gene")
    if new_feat:
        assert_files(output_gff, Path(data_dir / expected_gff))


@pytest.mark.parametrize(
    "in_gff, expected_gff",
    [
        param("ok_gene.gff", "ok_gene.gff", id="ok gene"),
        param("gene_ignored.gff", None, id="gene ignored"),
        param("mobile_te.gff", "mobile_te.gff", id="TE"),
        param("ok_protein_coding_gene.gff", "ok_gene.gff", id="ok protein_coding_gene"),
        param("ok_tr_ignored.gff", "ok_gene.gff", id="ok gene with ignored transcripts/subtranscripts"),
    ],
)
def test_simpler_gff3_feature(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: Optional[PathLike],
) -> None:
    """Test simplifying one gene (from a GFF3 file)."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
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
        param("mirna/mirna1.gff", "mirna/mirna1_simped.gff", does_not_raise(), id="miRNA split"),
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
    output_gff = tmp_dir / Path(in_gff).name
    with expectation:
        simp = GFFSimplifier()
        simp.simpler_gff3(input_gff)
        simp.records.to_gff(output_gff)
        assert_files(output_gff, data_dir / expected_gff)


@pytest.mark.parametrize(
    "in_gff, expected_gff, allow_cds",
    [
        param("ok_gene.gff", "ok_gene.gff", False, id="ok gene"),
        param("ok_gene.gff", "ok_gene.gff", True, id="ok gene, allow pseudo CDS"),
        param("pseudogene.gff", "pseudogene.gff", False, id="ok pseudogene"),
        param("pseudogene_cds.gff", "pseudogene_cds_removed.gff", False, id="pseudogene cds removed"),
        param("pseudogene_cds.gff", "pseudogene_cds.gff", True, id="pseudogene cds kept"),
    ],
)
def test_simpler_gff3_pseudogene(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    allow_cds: bool,
) -> None:
    """Test simplifying pseudogenes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / Path(in_gff).name
    simp = GFFSimplifier()
    simp.allow_pseudogene_with_cds = allow_cds
    simp.simpler_gff3(input_gff)
    simp.records.to_gff(output_gff)
    assert_files(output_gff, data_dir / expected_gff)


@pytest.mark.parametrize(
    "in_gff, expected_gff, skip_unrecognized, expectation",
    [
        param("bad_gene_type.gff", "", False, raises(GFFParserError), id="Unset skip unrecognized, fail"),
        param(
            "bad_gene_type.gff",
            "bad_gene_type_skipped.gff",
            True,
            does_not_raise(),
            id="True skip unrecognized, no fail",
        ),
        param("bad_gene_type.gff", "", False, raises(GFFParserError), id="bad type, fail"),
        param("ok_gene.gff", "ok_gene.gff", False, does_not_raise(), id="ok type, fail"),
    ],
)
def test_simpler_gff3_skip(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    skip_unrecognized: bool,
    expectation: ContextManager,
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    simp = GFFSimplifier(skip_unrecognized=skip_unrecognized)
    with expectation:
        simp.simpler_gff3(input_gff)
    if expected_gff:
        simp.records.to_gff(output_gff)
        assert_files(output_gff, data_dir / expected_gff)


@pytest.mark.parametrize(
    "genome_file, in_gff, expected_gff",
    [
        param(
            None,
            "genes_badnames.gff",
            "genes_badnames_noname.gff",
            id="Genes with bad names, no genome",
        ),
        param(
            "genome_no_brc4.json",
            "genes_badnames.gff",
            "genes_badnames_noname.gff",
            id="Genes with bad names, genome not BRC4",
        ),
        param(
            "genome_brc4.json",
            "genes_badnames.gff",
            "genes_badnames_brc4name.gff",
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
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / in_gff
    if genome_file is None:
        simp = GFFSimplifier()
    else:
        simp = GFFSimplifier(genome_path=data_dir / genome_file)
    simp.simpler_gff3(input_gff)
    simp.records.to_gff(output_gff)
    assert_files(output_gff, data_dir / expected_gff)


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param(
            "ok_gene.gff",
            "ok_gene.gff",
            does_not_raise(),
            id="Gene without miRNA",
        ),
        param(
            "mirna/mirna1.gff",
            "mirna/mirna1_simped.gff",
            does_not_raise(),
            id="Primary_transcript with 1 miRNA",
        ),
        param(
            "mirna/mirna2_pseudo.gff",
            "mirna/mirna2_pseudo_simped.gff",
            does_not_raise(),
            id="Pseudo primary_transcript",
        ),
        param(
            "mirna/mirna3_gene.gff",
            "mirna/mirna3_gene_simped.gff",
            does_not_raise(),
            id="Gene with Primary_transcript with 1 miRNA",
        ),
        param(
            "mirna/mirna4_unsupported.gff",
            "",
            raises(GFFParserError, match="Unknown subtype"),
            id="Gene with Primary_transcript with mRNA, not supported",
        ),
        param(
            "mirna/mirna5_unsupported.gff",
            "",
            raises(GFFParserError, match="too many sub_features"),
            id="Gene with Primary_transcript with 2 miRNAs, not supported",
        ),
    ],
)
def test_normalize_mirna(
    data_dir: Path,
    tmp_dir: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    expectation: ContextManager,
) -> None:
    """Test normalizing miRNA genes."""
    input_gff = data_dir / in_gff
    output_gff = tmp_dir / Path(in_gff).name
    with expectation:
        check_one_feature(input_gff, output_gff, "normalize_mirna")
        assert_files(output_gff, Path(data_dir / expected_gff))
