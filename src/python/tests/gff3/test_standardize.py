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
from typing import ContextManager, Optional

from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
from pytest import raises

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


def test_gene_to_cds(base_gene: SeqFeature):
    """Test the addition of intermediate transcripts."""

    cds = SeqFeature(base_gene.location, type="CDS")
    base_gene.sub_features.append(cds)

    gene = GFFStandard.gene_to_cds(base_gene)
    child = gene.sub_features[0]
    assert child.type == "mRNA"
    assert len(child.sub_features) == 2
