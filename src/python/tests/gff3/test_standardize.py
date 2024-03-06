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

from contextlib import nullcontext as does_not_raise
from typing import ContextManager, List, Optional

from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
from pytest import raises

from ensembl.io.genomio.gff3.exceptions import GFFParserError
from ensembl.io.genomio.gff3.standardize import GFFStandard


@pytest.fixture(name="base_gene")
def _base_gene() -> SeqFeature:
    gene_location = SimpleLocation(1, 100, 1)
    gene_source = "LOREM"
    gene = SeqFeature(gene_location, type="gene")
    gene.qualifiers["source"] = gene_source
    gene.sub_features = []
    return gene


def test_transcript_for_gene(base_gene: SeqFeature):
    """Test the creation of a transcript from a gene feature."""
    tr = GFFStandard.transcript_for_gene(base_gene)

    assert tr.location == base_gene.location
    assert tr.qualifiers["source"] == base_gene.qualifiers["source"]


@pytest.mark.dependency(name="build_transcript")
def test_build_transcript(
    base_gene: SeqFeature,
):
    """Test the creation of a transcript from a gene."""
    tr = GFFStandard.build_transcript(base_gene)
    assert tr.qualifiers["source"] == base_gene.qualifiers["source"]
    assert tr.sub_features == []


@pytest.mark.parametrize(
    "children, skip_non_cds, expected_mrna_children, expectation",
    [
        pytest.param(["CDS"], None, 2, does_not_raise(), id="One CDS"),
        pytest.param(["CDS", "CDS"], None, 4, does_not_raise(), id="Two CDS"),
        pytest.param(["exon"], None, 0, raises(GFFParserError), id="One exon"),
        pytest.param(["CDS", "exon"], True, 2, does_not_raise(), id="1 CDS + 1 exon, skip exon"),
    ],
)
@pytest.mark.dependency(name="gene_to_cds", depends=["build_transcript"])
def test_gene_to_cds(
    base_gene: SeqFeature,
    children: List[str],
    skip_non_cds: Optional[bool],
    expected_mrna_children: int,
    expectation: ContextManager,
):
    """Test the addition of intermediate transcripts."""

    for child in children:
        sub_feat = SeqFeature(base_gene.location, type=child)
        base_gene.sub_features.append(sub_feat)

    with expectation:
        if skip_non_cds is not None:
            gene = GFFStandard.gene_to_cds(base_gene, skip_non_cds=skip_non_cds)
        else:
            gene = GFFStandard.gene_to_cds(base_gene)

        # We should get only one mRNA and CDS + exons
        gene_child = gene.sub_features[0]
        assert gene_child.type == "mRNA"
        assert len(gene_child.sub_features) == expected_mrna_children