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
from typing import Any, ContextManager, Dict, List, Union

from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.exceptions import GFFParserError
from ensembl.io.genomio.gff3.standardize import (
    standardize_gene,
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
    """Generates features and structures for testing."""

    start = 1
    end = 1000
    strand = -1
    region = "LOREM"
    source = "Foo"

    def make(self, ftype: str, number: int = 1) -> List[SeqFeature]:
        """Create a defined number of features of a given type."""
        feats = []
        for _ in range(0, number):
            loc = SimpleLocation(self.start, self.end, self.strand)
            feat = SeqFeature(loc, type=ftype)
            feat.qualifiers["source"] = self.source
            feat.sub_features = []
            feats.append(feat)
        return feats

    def make_structure(self, children: List[Any]) -> List[SeqFeature]:
        """Return a SeqFeature children structure from the form:
        struct = ["mRNA"]
        struct = [{"mRNA": ["CDS", "exon"]}, "exon", "exon"]
        """
        output = []
        for child in children:
            if isinstance(child, str):
                output += self.make(child)
            elif isinstance(child, dict):
                feat_type = list(child.keys())[0]
                feat_children = list(child.values())[0]

                feat = self.make(feat_type)[0]
                feat.sub_features += self.make_structure(feat_children)
                output.append(feat)
            
        return output

    def get_sub_structure(self, feat: SeqFeature) -> Union[Dict, str]:
        """Create a children structure from a SeqFeature."""
        if feat.sub_features:
            feat_subs = []
            for sub in feat.sub_features:
                feat_subs.append(self.get_sub_structure(sub))
            return {feat.type: feat_subs}
        return feat.type


@pytest.mark.parametrize(
    "children, expected_children",
    [
        param(["gene"], [{"gene": ["transcript"]}], id="gene, no transcript: add"),
        param([{"gene": ["mRNA"]}], [{"gene": ["mRNA"]}], id="gene + mRNA"),
        param([{"gene": ["mRNA", "mRNA"]}], [{"gene": ["mRNA", "mRNA"]}], id="gene + 2 mRNA"),
        param(["ncRNA_gene"], ["ncRNA_gene"], id="1 ncRNA_gene, no transcript"),
        param([{"ncRNA_gene": ["transcript"]}], [{"ncRNA_gene": ["transcript"]}], id="1 ncRNA_gene + transcript"),
    ]
)
def test_add_transcript_to_naked_gene(children: List[Any], expected_children: List[Any]):
    """Test the creation of a transcript for a gene without one."""
    gen = FeatGenerator()
    genes = gen.make_structure(children)
    add_transcript_to_naked_gene(genes[0])
    assert gen.get_sub_structure(genes[0]) == expected_children[0]


@pytest.mark.parametrize(
    "children, expected_children, expectation",
    [
        param(["mRNA"], ["mRNA"], does_not_raise(), id="One mRNA, skip"),
        param(["CDS"], [{"mRNA": ["CDS", "exon"]}], does_not_raise(), id="One CDS, add mRNA"),
        param(["CDS", "CDS"], [{"mRNA": ["CDS", "exon", "CDS", "exon"]}], does_not_raise(), id="2 CDSs, add mRNA"),
        param(["exon"], ["exon"], does_not_raise(), id="One exon, skip"),
        param(["CDS", "exon"], ["CDS", "exon"], does_not_raise(), id="1 CDS + 1 exon, skip"),
    ],
)
def test_move_only_cdss_to_new_mrna(
    children: List[str],
    expected_children: Dict[str, int],
    expectation: ContextManager,
):
    """Test the creation of a new mRNA for CDSs under a gene."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)
    with expectation:
        move_only_cdss_to_new_mrna(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, expected_children, expectation",
    [
        param(["mRNA"], ["mRNA"], does_not_raise(), id="mRNA only, skip"),
        param(["CDS"], ["CDS"], does_not_raise(), id="1 CDS, skip"),
        param(["CDS", "exon"], ["CDS", "exon"], does_not_raise(), id="CDS + exon, skip"),
        param(["exon"], [{"mRNA": ["exon"]}], does_not_raise(), id="1 exon moved"),
        param(["exon", "exon"], [{"mRNA": ["exon", "exon"]}], does_not_raise(), id="2 exons moved"),
    ],
)
def test_move_only_exons_to_new_mrna(
    children: List[str],
    expected_children: Dict[str, int],
    expectation: ContextManager,
):
    """Test the creation of a new mRNA for exons under a gene."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)
    with expectation:
        move_only_exons_to_new_mrna(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, diff_exon, expected_children, expectation",
    [
        param(["mRNA"], False, ["mRNA"], does_not_raise(), id="mRNA only, skip"),
        param(["mRNA", "mRNA"], False, ["mRNA", "mRNA"], does_not_raise(), id="2 mRNA only, skip"),
        param(["CDS"], False, ["CDS"], does_not_raise(), id="One CDS, skip"),
        param(["CDS", "exon"], False, ["CDS", "exon"], does_not_raise(), id="CDS+exon, skip"),
        param(["mRNA", "CDS"], False, [{"mRNA": ["exon", "CDS"]}], does_not_raise(), id="1 mRNA + 1 CDS"),
        param(
            ["mRNA", "CDS", "extra"],
            False,
            ["extra", {"mRNA": ["exon", "CDS"]}],
            does_not_raise(),
            id="1 mRNA + 1 CDS + extra",
        ),
        param(
            ["mRNA", "CDS", "CDS"],
            False,
            [{"mRNA": ["exon", "exon", "CDS", "CDS"]}],
            does_not_raise(),
            id="1 mRNA + 2 CDS",
        ),
        param(["mRNA", "mRNA", "CDS"], False, [], raises(GFFParserError), id="2 mRNA only + CDS, fail"),
        param(
            ["mRNA", "mRNA", "CDS", "exon"],
            False,
            [],
            raises(GFFParserError),
            id="2 mRNA only + CDS + exon, fail",
        ),
        param([{"mRNA": ["CDS"]}, "CDS"], False, [], raises(GFFParserError), id="1 mRNA/CDS + 1 CDS, fail"),
        param(
            [{"mRNA": ["exon"]}, "CDS"],
            False,
            [{"mRNA": ["exon", "CDS"]}],
            does_not_raise(),
            id="1 mRNA/exon + 1 CDS same",
        ),
        param([{"mRNA": ["exon"]}, "CDS"], True, [], raises(GFFParserError), id="1 mRNA/exon + 1 CDS diff"),
        param(
            [{"mRNA": ["extra"]}, "CDS"],
            False,
            [{"mRNA": ["extra", "exon", "CDS"]}],
            does_not_raise(),
            id="1 mRNA/extra + 1 CDS",
        ),
        param(
            [{"mRNA": ["exon", "exon"]}, "CDS"],
            True,
            [],
            raises(GFFParserError),
            id="1 mRNA/2exon + 1 CDS same",
        ),
    ],
)
def test_move_cds_to_existing_mrna(
    children: List[str],
    diff_exon: bool,
    expected_children: Dict[str, int],
    expectation: ContextManager,
):
    """Test moving CDSs under a gene to under the mRNA."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)

    if diff_exon:
        for sub in gene.sub_features:
            if sub.type == "mRNA":
                for sub2 in sub.sub_features:
                    if sub2.type == "exon":
                        sub2.location += 10

    with expectation:
        move_cds_to_existing_mrna(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, has_id, expected_children, expectation",
    [
        param(["tRNA"], 0, ["tRNA"], does_not_raise(), id="tRNA only, skip"),
        param(["mRNA"], 0, ["mRNA"], does_not_raise(), id="mRNA only, skip"),
        param(["mRNA", "exon"], 0, ["mRNA", "exon"], does_not_raise(), id="mRNA and 1 exon without id"),
        param(["mRNA", "exon"], 1, ["mRNA"], does_not_raise(), id="mRNA and 1 exon with id"),
        param(
            [{"mRNA": ["exon"]}, "exon"],
            1,
            [{"mRNA": ["exon"]}],
            does_not_raise(),
            id="mRNA and 1 exon with id",
        ),
        param(
            [{"mRNA": ["exon"]}, "exon", "exon"],
            1,
            [],
            raises(GFFParserError),
            id="mRNA and 2 exons with partial id-",
        ),
        param(
            ["mRNA", "exon", "extra"],
            1,
            ["mRNA", "extra"],
            does_not_raise(),
            id="mRNA + extra + 1 exon with id",
        ),
    ],
)
def test_remove_extra_exons(
    children: List[Any], has_id: int, expected_children: List[Any], expectation: ContextManager
):
    """Test removing extra unneeded exons."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)

    if has_id:
        exon_num = 1
        for subfeat in gene.sub_features:
            if subfeat.type == "exon":
                subfeat.id = f"id-{exon_num}"
                exon_num += 1
            if exon_num > has_id:
                break

    with expectation:
        remove_extra_exons(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, expected_children, expectation",
    [
        param([{"mRNA": ["CDS", "exon"]}], [{"mRNA": ["CDS", "exon"]}], does_not_raise(), id="OK gene"),
        param(["mRNA", "CDS"], [{"mRNA": ["exon", "CDS"]}], does_not_raise(), id="mRNA + CDS, fixed"),
        param([{"mRNA": ["CDS", "exon"]}, "CDS"], [], raises(GFFParserError), id="Gene + extra CDS, fail"),
        param([{"mRNA": ["CDS", "exon"]}, "exon"], [], raises(GFFParserError), id="Gene + extra exon, fail"),
    ],
)
def test_standardize(children: List[Any], expected_children: List[Any], expectation: ContextManager) -> None:
    """Test the standardize() main function."""

    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)
    with expectation:
        standardize_gene(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}
