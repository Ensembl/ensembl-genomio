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
"""Unit testing of `ensembl.io.genomio.gff3.extract_annotation` module."""

from contextlib import nullcontext as does_not_raise
from typing import ContextManager, Dict, List, Optional

from Bio.SeqFeature import SeqFeature
import pytest
from pytest import raises, param

from ensembl.io.genomio.gff3.extract_annotation import (
    FunctionalAnnotations,
    MissingParentError,
    AnnotationError,
)


@pytest.mark.parametrize(
    "description, feature_id, output",
    [
        ("", [], False),
        ("", "PROTID12345", False),
        ("PROTID12345", ["PROTID12345"], False),
        ("ProtId12345", ["PROTID12345"], False),
        ("hypothetical PROTID12345 (ProtId12345)", ["PROTID12345"], False),
        ("hypothetical protein", [], False),
        ("hypothetical_protein", [], False),
        ("Hypothetical protein", [], False),
        ("hypothetical protein PROTID12345", ["PROTID12345"], False),
        ("hypothetical protein ProtId12345", ["PROTID12345"], False),
        ("hypothetical protein (ProtId12345)", ["PROTID12345"], False),
        ("hypothetical protein (fragment)", [], False),
        ("hypothetical protein, variant", [], False),
        ("hypothetical protein, variant 2", [], False),
        ("hypothetical protein - conserved", [], False),
        ("hypothetical protein, conserved", [], False),
        ("Hypothetical_protein_conserved", [], False),
        ("Hypothetical conserved protein", [], False),
        ("conserved hypothetical protein", [], False),
        ("conserved hypothetical protein, putative", [], False),
        ("conserved protein, unknown function", [], False),
        ("putative protein", [], False),
        ("putative_protein", [], False),
        ("hypothetical RNA", [], False),
        ("unspecified product", [], False),
        ("Unspecified product", [], False),
        ("conserved hypothetical transmembrane protein", [], True),
        ("unknown gene", [], False),
        ("unknown function", [], False),
        ("uncharacterized PROTID12345", ["PROTID12345"], False),
        ("LOW QUALITY PROTEIN: uncharacterized protein PROTID12345", ["PROTID12345"], False),
    ],
)
def test_product_is_informative(description: str, feature_id: Optional[List[str]], output: bool) -> None:
    """Tests the `FunctionalAnnotations.product_is_informative()` method."""
    assert FunctionalAnnotations.product_is_informative(description, feature_id) == output


@pytest.mark.parametrize(
    "seq_feat_type, feat_type, expected",
    [
        ("gene", "gene", does_not_raise()),
        ("mRNA", "transcript", does_not_raise()),
        ("CDS", "translation", does_not_raise()),
        ("transposable_element", "transposable_element", does_not_raise()),
        ("gene", "bad_type", raises(KeyError)),
    ],
)
@pytest.mark.dependency(name="add_feature")
def test_add_feature(seq_feat_type: str, feat_type: str, expected: ContextManager) -> None:
    """Tests the `FunctionaAnnotation.add_feature()` method with only one feature.

    Args:
        seq_feat_type: Type for the sequence feature to add.
        feat_type: Category type for that sequence feature.
        expected: What exception is expected to be raised, if any.

    """
    annot = FunctionalAnnotations()
    feature = SeqFeature(type=seq_feat_type, id="featA")
    with expected:
        annot.add_feature(feature, feat_type)
        assert annot.features[feat_type][feature.id]


@pytest.mark.parametrize(
    "parent_type, parent_id, child_id, expected",
    [
        ("gene", "geneA", "mrnA", does_not_raise()),
        ("bad_type", "geneA", "mrnA", raises(KeyError)),
        ("gene", "geneB", "mrnA", raises(MissingParentError)),
    ],
)
@pytest.mark.dependency(name="add_parent_link", depends=["add_feature"])
def test_add_parent_link(parent_type: str, parent_id: str, child_id: str, expected: ContextManager) -> None:
    """Tests the `FunctionaAnnotation.add_parent_link()` method.

    Add a parent feature, and then add a parent link.

    Args:
        parent_type: Type for the parent sequence feature.
        parent_id: ID for the parent sequence feature.
        child_id: ID for the child sequence feature.
        expected: What exception is expected to be raised, if any.

    """
    annot = FunctionalAnnotations()
    parent = SeqFeature(type="gene", id="geneA")
    annot.add_feature(parent, "gene")

    with expected:
        annot.add_parent_link(parent_type, parent_id, child_id)


@pytest.mark.parametrize(
    "in_parent_type, in_parent_id, in_child_id, out_parent_type, out_child_id, expected",
    [
        ("gene", "geneA", "mrnA", "gene", "mrnA", does_not_raise()),
        ("gene", "geneA", "mrnA", "bad_type", "mrnA", raises(KeyError)),
        ("gene", "geneA", "mrnA", "gene", "mrnB", raises(MissingParentError)),
    ],
)
@pytest.mark.dependency(name="get_parent", depends=["add_parent_link"])
def test_get_parent(
    in_parent_type: str,
    in_parent_id: str,
    in_child_id: str,
    out_parent_type: str,
    out_child_id: str,
    expected: ContextManager,
) -> None:
    """Tests the `FunctionaAnnotation.get_parent()` method.

    Args:
        in_parent_type: Type for the parent sequence feature.
        in_parent_id: ID for the parent sequence feature.
        in_child_id: ID for the child sequence feature.
        out_parent_type: Type for the parent stored in the functional annotation.
        out_child_id: ID for the child stored.
        expected: What exception is expected to be raised, if any.

    """
    annot = FunctionalAnnotations()
    parent = SeqFeature(type=in_parent_type, id=in_parent_id)
    annot.add_feature(parent, "gene")
    annot.add_feature(SeqFeature(type="mRNA", id=in_child_id), feat_type="transcript", parent_id=in_parent_id)

    with expected:
        out_parent = annot.get_parent(out_parent_type, out_child_id)
        assert out_parent == in_parent_id


@pytest.mark.parametrize(
    "child_type, child_id, out_parent_id, expected",
    [
        ("transcript", "mrna_A", "gene_A", does_not_raise()),
        pytest.param("bad_type", "mrna_A", "gene_A", raises(KeyError), id="Child type does not exist"),
        pytest.param("gene", "gene_A", None, raises(AnnotationError), id="Feature ID already loaded"),
        pytest.param(
            "gene", "gene_B", "gene_A", raises(AnnotationError), id="Cannot add a gene child of a gene"
        ),
    ],
)
@pytest.mark.dependency(name="add_feature_fail", depends=["add_feature", "get_parent"])
def test_add_feature_fail(
    child_type: str, child_id: str, out_parent_id: Optional[str], expected: ContextManager
) -> None:
    """Tests the `FunctionaAnnotation.add_feature()` method failures.

    Test the addition of a child feature after a parent has already been added.

    Args:
        child_type: Type for the child sequence feature.
        child_id: ID for the child sequence feature.
        out_parent_id: ID for the parent.
        expected: What exception is expected to be raised, if any.

    """
    annot = FunctionalAnnotations()
    parent = SeqFeature(type="gene", id="gene_A")
    child = SeqFeature(type="mRNA", id=child_id)
    annot.add_feature(parent, "gene")
    with expected:
        annot.add_feature(child, child_type, out_parent_id)


@pytest.mark.parametrize(
    "in_xrefs, genome, expected_xrefs",
    [
        param(None, None, [], id="No xref"),
        param([], None, [], id="Empty xref"),
        param(["DBname:Value"], None, [{"dbname": "DBname", "id": "Value"}], id="One xref"),
        param(
            ["DBname:Value:parts"],
            None,
            [{"dbname": "DBname", "id": "Value:parts"}],
            id="One xref with colon",
        ),
        param(["GO:XXX"], None, [], id="Ignore GO"),
        param(["GenBank:XXX"], None, [{"dbname": "GenBank", "id": "XXX"}], id="Genbank"),
        param(
            ["GenBank:XXX"],
            {"assembly": {"provider_name": "RefSeq"}},
            [{"dbname": "RefSeq", "id": "XXX"}],
            id="RefSeq",
        ),
        param(["GenBank:XXX"], {}, [{"dbname": "GenBank", "id": "XXX"}], id="Empty genome"),
        param(
            ["GenBank:XXX"],
            {"assembly": {}},
            [{"dbname": "GenBank", "id": "XXX"}],
            id="No provider in genome",
        ),
    ],
)
def test_get_xrefs(
    in_xrefs: Optional[List[str]], genome: Optional[Dict], expected_xrefs: List[Dict[str, str]]
) -> None:
    """Tests the `FunctionaAnnotation.get_xrefs()` method."""
    annot = FunctionalAnnotations()
    one_gene = SeqFeature(type="gene")
    if in_xrefs is not None:
        one_gene.qualifiers["Dbxref"] = in_xrefs
    if genome is not None:
        annot.genome = genome

    out_xrefs = annot.get_xrefs(one_gene)
    assert out_xrefs == expected_xrefs


@pytest.mark.parametrize(
    "feat_type, expected_number, expected",
    [
        ("gene", 1, does_not_raise()),
        ("transcript", 1, does_not_raise()),
        ("translation", 0, does_not_raise()),
        ("bad_type", 0, raises(KeyError)),
    ],
)
@pytest.mark.dependency(name="get_features", depends=["add_feature_fail"])
def test_get_features(feat_type: str, expected_number: int, expected: ContextManager) -> None:
    """Tests the `FunctionaAnnotation.get_features()` method.

    Load 2 features, then test the fetching of those features.

    Args:
        feat_type: Type for the features to fetch.
        expected: What exception is expected to be raised, if any.

    """
    annot = FunctionalAnnotations()
    one_gene = SeqFeature(type="gene", id="gene_A")
    one_transcript = SeqFeature(type="mRNA", id="mrna_A")
    annot.add_feature(one_gene, "gene")
    annot.add_feature(one_transcript, "transcript", parent_id=one_gene.id)

    with expected:
        out_feats = annot.get_features(feat_type)
        assert len(out_feats) == expected_number


@pytest.mark.parametrize(
    "gene_desc, transc_desc, transl_desc, out_gene_desc, out_transc_desc",
    [
        (None, None, None, None, None),
        ("Foobar", None, None, "Foobar", None),  # Only gene descriptions
        ("gene A", "transc B", "prod C", "gene A", "transc B"),  # All descriptions set
        (None, None, "Foobar", "Foobar", "Foobar"),  # Transfer from transl
        (None, "Foobar", None, "Foobar", "Foobar"),  # Transfer from transc
        (None, "Foobar", "Lorem", "Foobar", "Foobar"),  # Transfer from transc, transl also set
        ("Hypothetical gene", "Predicted function", "Foobar", "Foobar", "Foobar"),  # Non informative
        (None, None, "Unknown product", None, None),  # Non informative source
    ],
)
@pytest.mark.dependency(depends=["get_features"])
def test_transfer_descriptions(
    gene_desc: Optional[str],
    transc_desc: Optional[str],
    transl_desc: Optional[str],
    out_gene_desc: Optional[str],
    out_transc_desc: Optional[str],
) -> None:
    """Tests the `FunctionaAnnotation.transfer_descriptions()` method.

    Load 3 features (gene, transcript, translation) with or without a description for each one.

    Args:
        gene_desc: Description for the gene.
        transc_desc: Description for the transcript.
        transl_desc: Description for the translation.
        out_gene_desc: Excpected description for the gene after transfer.
        out_transc_desc: Excpected description for the transcript after transfer.

    """
    annot = FunctionalAnnotations()
    gene_name = "gene_A"
    transcript_name = "tran_A"
    one_gene = SeqFeature(type="gene", id=gene_name)
    if gene_desc:
        one_gene.qualifiers["description"] = [gene_desc]
    one_transcript = SeqFeature(type="mRNA", id=transcript_name)
    if transc_desc:
        one_transcript.qualifiers = {"product": [transc_desc]}
    one_translation = SeqFeature(type="CDS", id="cds_A")
    if transl_desc:
        one_translation.qualifiers = {"product": [transl_desc]}
    annot.add_feature(one_gene, "gene")
    annot.add_feature(one_transcript, "transcript", parent_id=one_gene.id)
    annot.add_feature(one_translation, "translation", parent_id=one_transcript.id)

    annot.transfer_descriptions()
    genes = annot.get_features("gene")
    transcs = annot.get_features("transcript")
    assert genes[gene_name].get("description") == out_gene_desc
    assert transcs[transcript_name].get("description") == out_transc_desc


@pytest.mark.dependency(depends=["add_feature"])
@pytest.mark.parametrize(
    "num_cds, cds_parts, expected_num_genes, expected_num_tr, expected_num_cds",
    [
        pytest.param(0, 0, 1, 1, 0, id="Store gene without CDS"),
        pytest.param(1, 1, 1, 1, 1, id="Store gene with 1 CDS in one part"),
        pytest.param(1, 2, 1, 1, 1, id="Store gene with 1 CDS in 2 parts"),
        pytest.param(2, 1, 1, 2, 2, id="Store gene with 2 CDS in 1 part each"),
        pytest.param(2, 2, 1, 2, 2, id="Store gene with 2 CDS in 2 part each"),
    ],
)
def test_store_gene(
    cds_parts: int, num_cds: int, expected_num_genes: int, expected_num_tr: int, expected_num_cds: int
) -> None:
    """Test store_gene given a gene Feature with a transcript and optional translation.

    Args:
        num_cds: Number of CDSs stored
        cds_parts: Number of parts of each CDS
        expected_num_genes: Number of genes stored as expected
        expected_num_tr: Number of transcripts stored as expected
        expected_num_cds: Number of CDSs stored as expected
    """
    annot = FunctionalAnnotations()
    gene_name = "gene_A"
    transcript_name = "tran_A"
    one_gene = SeqFeature(type="gene", id=gene_name)
    one_gene.sub_features = []

    # Add a translation (possibly in parts)
    if num_cds:
        for cds_number in range(1, num_cds + 1):
            transcript = SeqFeature(type="mRNA", id=f"tran_{cds_number}")
            transcript.sub_features = []
            exon = SeqFeature(type="exon", id=f"exon_{cds_number}")
            transcript.sub_features.append(exon)
            if cds_parts > 0:
                for _ in range(1, cds_parts + 1):
                    translation = SeqFeature(type="CDS", id=f"cds_{cds_number}")
                    transcript.sub_features.append(translation)
            one_gene.sub_features.append(transcript)
    else:
        one_transcript = SeqFeature(type="mRNA", id=transcript_name)
        one_transcript.sub_features = []
        one_exon = SeqFeature(type="exon", id="exon_A")
        one_transcript.sub_features.append(one_exon)
        one_gene.sub_features.append(one_transcript)

    annot.store_gene(one_gene)
    assert len(annot.features["gene"]) == expected_num_genes
    assert len(annot.features["transcript"]) == expected_num_tr
    assert len(annot.features["translation"]) == expected_num_cds
