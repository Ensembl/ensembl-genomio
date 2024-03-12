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
"""Unit testing of `ensembl.io.genomio.gff3.standardize` module."""

from collections import Counter
from contextlib import nullcontext as does_not_raise
from typing import ContextManager, Dict, List, Optional

from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.exceptions import GFFParserError
from ensembl.io.genomio.gff3.standardize import (
    # standardize_gene,
    add_transcript_to_naked_gene,
    move_only_cdss_to_new_mrna,
    move_only_exons_to_new_mrna,
    move_cds_to_existing_mrna,
    remove_extra_exons,
)


@pytest.fixture(name="base_gene")
def _base_gene() -> SeqFeature:
    gene_location = SimpleLocation(1, 100, 1)
    gene_source = "LOREM"
    gene = SeqFeature(gene_location, type="gene")
    gene.qualifiers["source"] = gene_source
    gene.sub_features = []
    return gene


class FeatGenerator:
    """Generates feature of a given type for testing."""

    start = 1
    end = 1000
    strand = -1
    region = "LOREM"
    source = "Foo"

    def make(self, ftype: str, number: int) -> List[SeqFeature]:
        """Create a defined bnumber of features of a given type."""
        feats = []
        for _ in range(0, number):
            loc = SimpleLocation(self.start, self.end, self.strand)
            feat = SeqFeature(loc, type=ftype)
            feat.qualifiers["source"] = self.source
            feat.sub_features = []
            feats.append(feat)
        return feats

    def append(self, feat: SeqFeature, ftype: str, number: int) -> None:
        """Create a defined bnumber of features of a given type and append them to the gene."""
        subs = self.make(ftype, number)
        feat.sub_features += subs


@pytest.mark.parametrize(
    "gene_type, ntr_before, ntr_after",
    [
        param("ncRNA_gene", 1, 1, id="ncRNA_gene with 1 transcript"),
        param("ncRNA_gene", 0, 0, id="ncRNA_gene with no transcript"),
        param("gene", 1, 1, id="Gene with 1 transcript"),
        param("gene", 2, 2, id="Gene with 2 transcripts"),
        param("gene", 0, 1, id="Gene with no transcript"),
    ],
)
def test_add_transcript_to_naked_gene(gene_type: str, ntr_before: int, ntr_after: int):
    """Test the creation of a transcript from a gene feature."""
    gen = FeatGenerator()
    gene = gen.make(gene_type, 1)[0]
    if ntr_before > 0:
        gen.append(gene, "mRNA", ntr_before)

    add_transcript_to_naked_gene(gene)
    assert len(gene.sub_features) == ntr_after


@pytest.mark.parametrize(
    "children, expected_children, expected_mrna_children, expectation",
    [
        pytest.param(["mRNA"], {"mRNA": 1}, {}, does_not_raise(), id="One mRNA, skip"),
        pytest.param(["CDS"], {"mRNA": 1}, {"CDS": 1, "exon": 1}, does_not_raise(), id="One CDS, add mRNA"),
        pytest.param(
            ["CDS", "CDS"], {"mRNA": 1}, {"CDS": 2, "exon": 2}, does_not_raise(), id="Two CDS, add mRNA"
        ),
        pytest.param(["exon"], {"exon": 1}, {}, does_not_raise(), id="One exon, skip"),
        pytest.param(["CDS", "exon"], {"CDS": 1, "exon": 1}, {}, does_not_raise(), id="1 CDS + 1 exon, skip"),
    ],
)
def test_move_only_cdss_to_new_mrna(
    children: List[str],
    expected_children: Dict[str, int],
    expected_mrna_children: Dict[str, int],
    expectation: ContextManager,
):
    """Test the addition of intermediate transcripts."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]

    for child_type in children:
        gen.append(gene, child_type, 1)

    with expectation:
        move_only_cdss_to_new_mrna(gene)

        fcounter = dict(Counter([feat.type for feat in gene.sub_features]))
        assert fcounter == expected_children

        gene_child = gene.sub_features[0]
        if gene_child.type == "mRNA":
            fcounter_mrna = dict(Counter([feat.type for feat in gene_child.sub_features]))
            assert fcounter_mrna == expected_mrna_children


@pytest.mark.parametrize(
    "children, expected_children, expected_mrna_children, expectation",
    [
        pytest.param(["mRNA"], {"mRNA": 1}, {}, does_not_raise(), id="mRNA only, skip"),
        pytest.param(["CDS"], {"CDS": 1}, {}, does_not_raise(), id="One CDS, skip"),
        pytest.param(["CDS", "exon"], {"CDS": 1, "exon": 1}, {}, does_not_raise(), id="One CDS, skip"),
        pytest.param(["exon"], {"mRNA": 1}, {"exon": 1}, does_not_raise(), id="One exon"),
        pytest.param(["exon", "exon"], {"mRNA": 1}, {"exon": 2}, does_not_raise(), id="Two exons"),
    ],
)
def test_move_only_exons_to_new_mrna(
    children: List[str],
    expected_children: Dict[str, int],
    expected_mrna_children: Dict[str, int],
    expectation: ContextManager,
):
    """Test the addition of intermediate transcripts."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]

    for child_type in children:
        gen.append(gene, child_type, 1)

    with expectation:
        move_only_exons_to_new_mrna(gene)

        fcounter = dict(Counter([feat.type for feat in gene.sub_features]))
        assert fcounter == expected_children

        gene_child = gene.sub_features[0]
        if gene_child.type == "mRNA":
            fcounter_mrna = dict(Counter([feat.type for feat in gene_child.sub_features]))
            assert fcounter_mrna == expected_mrna_children


@pytest.mark.parametrize(
    "children, expected_children, expected_mrna_children, expectation",
    [
        pytest.param(["mRNA"], {"mRNA": 1}, {}, does_not_raise(), id="mRNA only, skip"),
        pytest.param(["mRNA", "mRNA"], {"mRNA": 2}, {}, does_not_raise(), id="2 mRNA only, skip"),
        pytest.param(["CDS"], {"CDS": 1}, {}, does_not_raise(), id="One CDS, skip"),
        pytest.param(["CDS", "exon"], {"CDS": 1, "exon": 1}, {}, does_not_raise(), id="CDS+exon, skip"),
        pytest.param(
            ["mRNA", "CDS"], {"mRNA": 1}, {"CDS": 1, "exon": 1}, does_not_raise(), id="1 mRNA + 1 CDS"
        ),
        pytest.param(
            ["mRNA", "CDS", "CDS"], {"mRNA": 1}, {"CDS": 2, "exon": 2}, does_not_raise(), id="1 mRNA + 2 CDS"
        ),
        pytest.param(["mRNA", "mRNA", "CDS"], {}, {}, raises(GFFParserError), id="2 mRNA only + CDS, fail"),
        pytest.param(
            ["mRNA", "mRNA", "CDS", "exon"],
            {"mRNA": 2, "CDS": 1, "exon": 1},
            {},
            raises(GFFParserError),
            id="2 mRNA only + CDS + exon, fail",
        ),
    ],
)
def test_move_cds_to_existing_mrna(
    children: List[str],
    expected_children: Dict[str, int],
    expected_mrna_children: Dict[str, int],
    expectation: ContextManager,
):
    """Test the addition of intermediate transcripts."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]

    for child_type in children:
        gen.append(gene, child_type, 1)

    with expectation:
        move_cds_to_existing_mrna(gene)

        fcounter = dict(Counter([feat.type for feat in gene.sub_features]))
        assert fcounter == expected_children

        gene_child = gene.sub_features[0]
        if gene_child.type == "mRNA":
            fcounter_mrna = dict(Counter([feat.type for feat in gene_child.sub_features]))
            assert fcounter_mrna == expected_mrna_children


def test_remove_extra_exons():
    """Test the addition of intermediate transcripts."""
    # TODO
