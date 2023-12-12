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
"""Unit testing of `ensembl.io.genomio.gff3.extract_annotation` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_extract_annotation.py

"""

from typing import Optional
import pytest

from Bio.SeqFeature import SeqFeature

from ensembl.io.genomio.gff3.extract_annotation import FunctionalAnnotations


class TestFunctionalAnnotations:
    """Tests for the `FunctionalAnnotations` class."""

    @pytest.mark.parametrize(
        "feature, feat_type, parent_id",
        [
            (SeqFeature(type="gene", id="geneA"), "gene", None),
            (SeqFeature(type="mRNA", id="mrnaA"), "transcript", None),
            (SeqFeature(type="CDS", id="cdsA"), "translation", None),
            (SeqFeature(type="transposable_element", id="teA"), "transposable_element", None),
        ]
    )
    def test_add_feature(self, feature: SeqFeature, feat_type: str, parent_id: Optional[str]):
        annot = FunctionalAnnotations()
        annot.add_feature(feature, feat_type, parent_id)
        assert annot.features[feat_type][feature.id]

    def test_add_parent(self):
        parent = SeqFeature(type="gene", id="geneA")
        child = SeqFeature(type="transcript", id="mrnaA")
        annot = FunctionalAnnotations()
        annot.add_feature(parent, "gene")
        annot.add_parent("gene", parent.id, child.id)

        assert annot.parents["gene"][child.id] == parent.id

    def test_add_feature_parents(self):
        parent = SeqFeature(type="gene", id="geneA")
        child = SeqFeature(type="transcript", id="mrnaA")
        annot = FunctionalAnnotations()
        annot.add_feature(parent, "gene")
        annot.add_feature(child, "transcript", parent_id=parent.id)

        assert annot.parents["gene"][child.id] == parent.id

    def test_add_feature_get_parents(self):
        parent = SeqFeature(type="gene", id="geneA")
        child = SeqFeature(type="transcript", id="mrnaA")
        annot = FunctionalAnnotations()
        annot.add_feature(parent, "gene")
        annot.add_feature(child, "transcript", parent_id=parent.id)

        parent_id = annot.get_parent("gene", child.id)

        assert parent_id == parent.id

    @pytest.mark.parametrize(
        "description, feature_id, output",
        [
            ("", None, False),
            ("", "PROTID12345", False),
            ("PROTID12345", "PROTID12345", False),
            ("ProtId12345", "PROTID12345", False),
            ("hypothetical PROTID12345 (ProtId12345)", "PROTID12345", False),
            ("hypothetical protein", None, False),
            ("hypothetical_protein", None, False),
            ("Hypothetical protein", None, False),
            ("hypothetical protein PROTID12345", "PROTID12345", False),
            ("hypothetical protein ProtId12345", "PROTID12345", False),
            ("hypothetical protein (ProtId12345)", "PROTID12345", False),
            ("hypothetical protein (fragment)", None, False),
            ("hypothetical protein, variant", None, False),
            ("hypothetical protein, variant 2", None, False),
            ("hypothetical protein - conserved", None, False),
            ("hypothetical protein, conserved", None, False),
            ("Hypothetical_protein_conserved", None, False),
            ("Hypothetical conserved protein", None, False),
            ("conserved hypothetical protein", None, False),
            ("conserved hypothetical protein, putative", None, False),
            ("conserved protein, unknown function", None, False),
            ("putative protein", None, False),
            ("putative_protein", None, False),
            ("hypothetical RNA", None, False),
            ("unspecified product", None, False),
            ("Unspecified product", None, False),
            ("conserved hypothetical transmembrane protein", None, True),
        ],
    )
    def test_product_is_informative(self, description: str, feature_id: Optional[str], output: bool) -> None:
        """Tests the `FunctionalAnnotations.product_is_informative()` method.

        Args:
            description: Product description.
            feature_id: Feature ID that might be in the description.
            output: Expected output, i.e. True if description is informative, False otherwise.

        """
        assert FunctionalAnnotations.product_is_informative(description, feature_id) == output
