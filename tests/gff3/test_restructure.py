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
from typing import Any, ContextManager, Dict, List, Union

from Bio.SeqFeature import SimpleLocation
import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.exceptions import GFFParserError
from ensembl.io.genomio.gff3 import restructure
from ensembl.io.genomio.gff3.features import GFFSeqFeature


class FeatGenerator:
    """Generates features and structures for testing."""

    start = 1
    end = 1000
    strand = -1
    region = "LOREM"
    source = "Foo"

    def make(self, ftype: str, number: int = 1) -> List[GFFSeqFeature]:
        """Returns a list with a defined number of features of a given type."""
        feats = []
        for _ in range(0, number):
            loc = SimpleLocation(self.start, self.end, self.strand)
            feat = GFFSeqFeature(loc, type=ftype)
            feat.qualifiers["source"] = self.source
            feats.append(feat)
        return feats

    def make_structure(self, children: List[Any]) -> List[GFFSeqFeature]:
        """Returns a list of SeqFeature children structure from the form:
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

    def get_sub_structure(self, feat: GFFSeqFeature) -> Union[Dict, str]:
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
        param(
            [{"ncRNA_gene": ["transcript"]}], [{"ncRNA_gene": ["transcript"]}], id="1 ncRNA_gene + transcript"
        ),
    ],
)
def test_add_transcript_to_naked_gene(children: List[Any], expected_children: List[Any]) -> None:
    """Test the creation of a transcript for a gene without one."""
    gen = FeatGenerator()
    genes = gen.make_structure(children)
    restructure.add_transcript_to_naked_gene(genes[0])
    assert gen.get_sub_structure(genes[0]) == expected_children[0]


@pytest.mark.parametrize(
    "children, expected_children",
    [
        param(["mRNA"], ["mRNA"], id="One mRNA, skip"),
        param(["CDS"], [{"mRNA": ["CDS", "exon"]}], id="One CDS, add mRNA"),
        param(
            ["CDS", "CDS"],
            [{"mRNA": ["CDS", "exon", "CDS", "exon"]}],
            id="2 CDSs, add mRNA",
        ),
        param(["exon"], ["exon"], id="One exon, skip"),
        param(["CDS", "exon"], ["CDS", "exon"], id="1 CDS + 1 exon, skip"),
    ],
)
def test_move_only_cdss_to_new_mrna(
    children: List[str],
    expected_children: Dict[str, int],
) -> None:
    """Test the creation of a new mRNA for CDSs under a gene."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)
    restructure.move_only_cdss_to_new_mrna(gene)
    assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, expected_children",
    [
        param(["mRNA"], ["mRNA"], id="mRNA only, skip"),
        param(["CDS"], ["CDS"], id="1 CDS, skip"),
        param(["CDS", "exon"], ["CDS", "exon"], id="CDS + exon, skip"),
        param(["exon"], [{"mRNA": ["exon"]}], id="1 exon moved"),
        param(["exon", "exon"], [{"mRNA": ["exon", "exon"]}], id="2 exons moved"),
    ],
)
def test_move_only_exons_to_new_mrna(
    children: List[str],
    expected_children: Dict[str, int],
) -> None:
    """Test the creation of a new mRNA for exons under a gene."""
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)
    restructure.move_only_exons_to_new_mrna(gene)
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
        param(
            ["mRNA", "mRNA", "CDS"],
            False,
            [],
            raises(GFFParserError, match="contains several mRNAs"),
            id="2 mRNA only + CDS, fail",
        ),
        param(
            ["mRNA", "mRNA", "CDS", "exon"],
            False,
            [],
            raises(GFFParserError, match="contains several mRNAs"),
            id="2 mRNA only + CDS + exon, fail",
        ),
        param(
            [{"mRNA": ["CDS"]}, "CDS"],
            False,
            [],
            raises(GFFParserError, match="children in both"),
            id="1 mRNA/CDS + 1 CDS, fail",
        ),
        param(
            [{"mRNA": ["exon"]}, "CDS"],
            False,
            [{"mRNA": ["exon", "CDS"]}],
            does_not_raise(),
            id="1 mRNA/exon + 1 CDS same",
        ),
        param(
            [{"mRNA": ["exon"]}, "CDS"],
            True,
            [],
            raises(GFFParserError, match="do not match"),
            id="1 mRNA/exon + 1 CDS diff",
        ),
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
            raises(GFFParserError, match="different count"),
            id="1 mRNA/2exon + 1 CDS same",
        ),
    ],
)
def test_move_cds_to_existing_mrna(
    children: List[str],
    diff_exon: bool,
    expected_children: Dict[str, int],
    expectation: ContextManager,
) -> None:
    """Test moving CDSs under a gene to under the mRNA.

    Args:
        diff_exons: use exons with different coordinates than the CDSs.
    """
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
        restructure.move_cds_to_existing_mrna(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, has_id, expected_children, expectation",
    [
        param(["tRNA"], 0, ["tRNA"], does_not_raise(), id="tRNA only, skip"),
        param(["mRNA"], 0, ["mRNA"], does_not_raise(), id="mRNA only, skip"),
        param(
            ["mRNA", "exon"],
            0,
            ["mRNA", "exon"],
            raises(GFFParserError, match="not all start"),
            id="mRNA and 1 exon without id",
        ),
        param(["mRNA", "exon"], 1, ["mRNA"], does_not_raise(), id="mRNA and 1 exon with id"),
        param(
            [{"mRNA": ["exon"]}, "exon"],
            1,
            [{"mRNA": ["exon"]}],
            does_not_raise(),
            id="mRNA with exon and 1 exon with id",
        ),
        param(
            [{"mRNA": ["exon"]}, "exon", "exon"],
            1,
            [],
            raises(GFFParserError, match="not all start"),
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
) -> None:
    """Test removing extra unneeded exons.

    Args:
        has_id: add an ID starting with 'id-' for this number of exons (if any).
    """
    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)

    if has_id:
        exon_num = 1
        for subfeature in gene.sub_features:
            if subfeature.type == "exon":
                subfeature.id = f"id-{exon_num}"
                exon_num += 1
            if exon_num > has_id:
                break

    with expectation:
        restructure.remove_extra_exons(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, expected_children, expectation",
    [
        param([{"mRNA": ["CDS", "exon"]}], [{"mRNA": ["CDS", "exon"]}], does_not_raise(), id="OK gene"),
        param(["mRNA", "CDS"], [{"mRNA": ["exon", "CDS"]}], does_not_raise(), id="mRNA + CDS, fixed"),
        param(
            [{"mRNA": ["CDS", "exon"]}, "CDS"],
            [],
            raises(GFFParserError, match="in both"),
            id="Gene + extra CDS, fail",
        ),
        param(
            [{"mRNA": ["CDS", "exon"]}, "exon"],
            [],
            raises(GFFParserError, match="not all start"),
            id="Gene + extra exon, fail",
        ),
    ],
)
def test_restructure_gene(
    children: List[Any], expected_children: List[Any], expectation: ContextManager
) -> None:
    """Test the `restructure_gene()` main function."""

    gen = FeatGenerator()
    gene = gen.make("gene", 1)[0]
    gene.sub_features += gen.make_structure(children)
    with expectation:
        restructure.restructure_gene(gene)
        assert gen.get_sub_structure(gene) == {"gene": expected_children}


@pytest.mark.parametrize(
    "children, expected_children",
    [
        param("gene", "gene", id="gene"),
        param("pseudogene", "pseudogene", id="pseudogene"),
        param({"pseudogene": ["mRNA"]}, {"pseudogene": ["mRNA"]}, id="pseudogene mRNA"),
        param(
            {"pseudogene": [{"mRNA": ["CDS", "CDS"]}]}, {"pseudogene": ["mRNA"]}, id="pseudogene mRNA CDSs"
        ),
        param(
            {"pseudogene": [{"mRNA": ["CDS", "exon"]}]},
            {"pseudogene": [{"mRNA": ["exon"]}]},
            id="pseudogene mRNA CDSs, exons",
        ),
        param({"pseudogene": ["CDS", "CDS"]}, "pseudogene", id="pseudogene CDSs"),
    ],
)
def test_remove_cds_from_pseudogene(children: List[Any], expected_children: List[Any]) -> None:
    """Test CDS removal from pseudogene."""
    gen = FeatGenerator()
    gene = gen.make_structure([children])[0]
    restructure.remove_cds_from_pseudogene(gene)
    assert gen.get_sub_structure(gene) == expected_children
