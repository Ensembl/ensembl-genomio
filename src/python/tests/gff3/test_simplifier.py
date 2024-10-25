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
"""Unit testing of `ensembl.io.genomio.gff3.simplifier` module."""
# pylint: disable=too-many-positional-arguments

from contextlib import nullcontext as does_not_raise
from os import PathLike
from pathlib import Path
from typing import Callable, ContextManager, Dict, Optional

import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.exceptions import GeneSegmentError, GFFParserError
from ensembl.io.genomio.gff3.simplifier import GFFSimplifier
from ensembl.io.genomio.gff3.exceptions import IgnoredFeatureError, UnsupportedFeatureError
from ensembl.io.genomio.gff3.features import GFFSeqFeature
from ensembl.io.genomio.utils import print_json


@pytest.mark.parametrize(
    "genome_meta, expected_provider_name",
    [
        param({}, "GenBank", id="No metadata"),
        param(
            {"assembly": {"provider_name": "LOREM", "accession": "GCA000"}},
            "LOREM",
            id="Explicit provider name",
        ),
        param({"assembly": {"accession": "GCA00000"}}, "GenBank", id="Genbank from accession"),
        param({"assembly": {"accession": "GCF00000"}}, "RefSeq", id="RefSeq from accession"),
        param(
            {"assembly": {"provider_name": "LOREM", "accession": "GCA00000"}},
            "LOREM",
            id="Explicit provider_name, GCA accession",
        ),
        param(
            {"assembly": {"provider_name": "LOREM", "accession": "GCF00000"}},
            "LOREM",
            id="Explicit provider_name, GCF accession",
        ),
    ],
)
def test_get_provider_name(tmp_path: Path, genome_meta: Dict, expected_provider_name: str) -> None:
    """Tests `GFFSimplifier.get_provider_name().`"""
    # Write metadata file
    meta_path = tmp_path / "meta.json"
    print_json(meta_path, genome_meta)
    simp = GFFSimplifier(meta_path)
    assert simp.get_provider_name() == expected_provider_name


@pytest.mark.parametrize(
    "genome_meta, expected_provider_name",
    [
        param({}, "GenBank", id="No metadata"),
        param(
            {"assembly": {"provider_name": "LOREM", "accession": "GCA000"}},
            "LOREM",
            id="Explicit provider name",
        ),
    ],
)
def test_init_provider_name(tmp_path: Path, genome_meta: Dict, expected_provider_name: str) -> None:
    """Tests `GFFSimplifier.__init__` to set the `provider_name` to its `FunctionalAnnotations` attrib."""
    # Write metadata file
    meta_path = tmp_path / "meta.json"
    print_json(meta_path, genome_meta)
    simp = GFFSimplifier(meta_path)

    assert simp.annotations.provider_name == expected_provider_name


def check_one_feature(input_gff: PathLike, output_gff: PathLike, check_function: str) -> None:
    """Load 1 feature from a GFF, apply a function, then write it back to a GFF."""
    simp = GFFSimplifier()
    simp.records.from_gff(input_gff)
    # Get the only feature
    feat = simp.records[0].features[0]
    # Apply the named function
    check_method = getattr(simp, check_function)
    new_feat = check_method(feat)
    # Put it back
    if isinstance(new_feat, list):
        simp.records[0].features = new_feat
    else:
        simp.records[0].features = [new_feat]
    simp.records.to_gff(output_gff)


@pytest.mark.parametrize(
    "in_gff, expected_gff",
    [
        param("ok_gene.gff", "ok_gene.gff", id="ok gene"),
        param("lone/transcript.gff", "lone/transcript_simped.gff", id="lone transcript"),
        param("lone/trna.gff", "lone/trna_simped.gff", id="lone tRNA"),
        param("lone/rrna.gff", "lone/rrna_simped.gff", id="lone rRNA"),
        param("lone/mrna.gff", "lone/mrna_simped.gff", id="lone mRNA"),
        param("lone/mrna_pseudo.gff", "lone/mrna_pseudo_simped.gff", id="lone pseudo mRNA"),
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
    check_one_feature(input_gff, output_gff, "create_gene_for_lone_transcript")
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
    tmp_path: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
) -> None:
    """Test gene create gene for lone CDS."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / Path(in_gff).name
    check_one_feature(input_gff, output_gff, "create_gene_for_lone_cds")
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
    feat = GFFSeqFeature(None, type=in_type)
    feat.qualifiers = {"source": "LOREM"}
    if in_mobile_type is not None:
        feat.qualifiers["mobile_element_type"] = [in_mobile_type]
    if in_product is not None:
        feat.qualifiers["product"] = [in_product]
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
    feat = GFFSeqFeature(None, type="non_gene_name")
    with raises(NotImplementedError):
        simp.normalize_non_gene(feat)


@pytest.mark.parametrize(
    "in_type, tr_name, out_type, expectation",
    [
        param("mRNA", "", "mRNA", does_not_raise(), id="mRNA no change"),
        param("C_gene_segment", "", "C_gene_segment", raises(GeneSegmentError), id="no standard name"),
        param("C_gene_segment", "immunoglobulin", "IG_C_gene", does_not_raise(), id="C immunoglobulin"),
        param("C_gene_segment", "ig", "IG_C_gene", does_not_raise(), id="C ig"),
        param("V_gene_segment", "t-cell", "TR_V_gene", does_not_raise(), id="V t-cell"),
        param("V_gene_segment", "T_cell", "TR_V_gene", does_not_raise(), id="V T_cell"),
        param("V_gene_segment", "Lorem Ipsum", "", raises(GeneSegmentError), id="V T_cell"),
    ],
)
def test_format_gene_segments(
    in_type: str,
    tr_name: str,
    out_type: str,
    expectation: ContextManager,
) -> None:
    """Test `format_gene_segments` without a CDS."""
    simp = GFFSimplifier()
    feat = GFFSeqFeature(None, type=in_type)
    if tr_name:
        feat.qualifiers["standard_name"] = [tr_name]
    with expectation:
        new_feat = simp.format_gene_segments(feat)
        assert new_feat.type == out_type


@pytest.mark.parametrize(
    "has_cds, cds_name, out_type, expectation",
    [
        param(False, "", "", raises(GeneSegmentError), id="No CDS"),
        param(True, "", "", raises(GeneSegmentError), id="CDS no info"),
        param(True, "ig", "IG_C_gene", does_not_raise(), id="C ig"),
        param(True, "t-cell", "TR_C_gene", does_not_raise(), id="C t-cell"),
    ],
)
def test_format_gene_segments_cds(
    has_cds: bool,
    cds_name: str,
    out_type: str,
    expectation: ContextManager,
) -> None:
    """Test `format_gene_segments` with a CDS (and no info on the transcript)."""
    simp = GFFSimplifier()
    feat = GFFSeqFeature(None, type="C_gene_segment")
    if has_cds:
        cds = GFFSeqFeature(None, type="CDS")
        cds.qualifiers["product"] = [cds_name]
        feat.sub_features = [cds]
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
    tmp_path: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
) -> None:
    """Test clean gene."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / Path(in_gff).name
    check_one_feature(input_gff, output_gff, "clean_gene")
    assert_files(output_gff, Path(data_dir / expected_gff))


@pytest.mark.parametrize(
    "in_gff, expected_gff, expectation",
    [
        param("ok_gene.gff", "ok_gene.gff", does_not_raise(), id="ok gene"),
        param("gene_ignored.gff", None, raises(IgnoredFeatureError), id="gene ignored"),
        param("gene_unsupported.gff", None, raises(UnsupportedFeatureError), id="gene unsupported"),
        param("mobile_te.gff", "mobile_te.gff", does_not_raise(), id="TE"),
        param("ok_protein_coding_gene.gff", "ok_gene.gff", does_not_raise(), id="ok protein_coding_gene"),
        param(
            "ok_tr_ignored.gff",
            "ok_gene.gff",
            does_not_raise(),
            id="ok gene with ignored transcripts/subtranscripts",
        ),
    ],
)
def test_simpler_gff3_feature(
    data_dir: Path,
    tmp_path: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: Optional[PathLike],
    expectation: ContextManager,
) -> None:
    """Test simplifying one gene (from a GFF3 file)."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / in_gff
    with expectation:
        check_one_feature(input_gff, output_gff, "simpler_gff3_feature")
    if expected_gff is not None:
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
    tmp_path: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    expectation: ContextManager,
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / Path(in_gff).name
    simp = GFFSimplifier()
    with expectation:
        simp.simpler_gff3(input_gff)
    if expected_gff:
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
    tmp_path: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    allow_cds: bool,
) -> None:
    """Test simplifying pseudogenes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / Path(in_gff).name
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
    tmp_path: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    skip_unrecognized: bool,
    expectation: ContextManager,
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / in_gff
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
    tmp_path: Path,
    assert_files: Callable,
    genome_file: Optional[PathLike],
    in_gff: PathLike,
    expected_gff: PathLike,
) -> None:
    """Test simplifying genes from GFF3 files."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / in_gff
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
        param("ok_gene.gff", "ok_gene.gff", does_not_raise(), id="normal gene"),
        param(
            "mirna/gene.gff",
            "mirna/gene_simped.gff",
            does_not_raise(),
            id="gene + primary_transcript + miRNA",
        ),
        param(
            "mirna/pseudogene.gff",
            "mirna/pseudogene_simped.gff",
            does_not_raise(),
            id="gene + primary_transcript - miRNA",
        ),
        param(
            "mirna/nogene.gff",
            "mirna/nogene_simped.gff",
            does_not_raise(),
            id="primary_transcript + miRNA",
        ),
        param(
            "mirna/pseudo_nogene.gff",
            "mirna/pseudo_nogene_simped.gff",
            does_not_raise(),
            id="primary_transcript - miRNA",
        ),
        param(
            "mirna/unsupported_tr.gff",
            "",
            raises(GFFParserError, match="Unknown subtype"),
            id="gene + primary_transcript + mRNA, not supported",
        ),
        param(
            "mirna/two_primary.gff",
            "",
            raises(GFFParserError, match="too many sub_features"),
            id="gene + 2x primary_transcript, not supported",
        ),
    ],
)
def test_simpler_gff3_mirna(
    data_dir: Path,
    tmp_path: Path,
    assert_files: Callable,
    in_gff: PathLike,
    expected_gff: PathLike,
    expectation: ContextManager,
) -> None:
    """Test normalizing miRNA genes."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / Path(in_gff).name
    simp = GFFSimplifier()
    with expectation:
        simp.simpler_gff3(input_gff)
    if expected_gff:
        simp.records.to_gff(output_gff)
        assert_files(output_gff, data_dir / expected_gff)
