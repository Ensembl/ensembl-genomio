#!/usr/bin/env python
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
"""Unit testing of :mod:`ensembl.io.genomio.schemas` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_schemas.py

"""

import pytest

from ensembl.io.genomio.gff3.functional_annotation import FunctionalAnnotations as fa


class TestMergeGFF3:
    """Tests for the integrity module."""

    @pytest.mark.parametrize(
        "description",
        [
            "",
            "hypothetical protein",
            "hypothetical_protein",
            "Hypothetical protein",
            "hypothetical protein (fragment)",
            "hypothetical protein, variant",
            "hypothetical protein, variant 2",
            "hypothetical protein - conserved",
            "hypothetical protein, conserved",
            "Hypothetical_protein_conserved",
            "Hypothetical conserved protein",
            "conserved hypothetical protein",
            "conserved hypothetical protein, putative",
            "conserved protein, unknown function",
            "putative protein",
            "putative_protein",
            "hypothetical RNA",
            "unspecified product",
            "Unspecified product",
        ],
    )
    def test_invalid_description(self, description: str) -> None:
        """Tests `functional_annotation.product_is_informative` method.

        Args:
            description: Description string to check.

        """
        assert not fa.product_is_informative(description)

    @pytest.mark.parametrize(
        "description",
        [
            "conserved hypothetical transmembrane protein",
        ],
    )
    def test_valid_description(self, description: str) -> None:
        """Tests `functional_annotation.product_is_informative` method with valid descriptions.

        Args:
            description: Description string to check.

        """
        assert fa.product_is_informative(description)

    @pytest.mark.parametrize(
        "description, feat_id",
        [
            ("", "PROTID12345"),
            ("PROTID12345", "PROTID12345"),
            ("ProtId12345", "PROTID12345"),
            ("hypothetical protein PROTID12345", "PROTID12345"),
            ("hypothetical protein ProtId12345", "PROTID12345"),
            ("hypothetical protein (ProtId12345)", "PROTID12345"),
            ("hypothetical PROTID12345 (ProtId12345)", "PROTID12345"),
        ],
    )
    def test_invalid_description_with_id(self, description: str, feat_id: str) -> None:
        """Tests `functional_annotation.product_is_informative` method with an ID in the description.

        Args:
            description: Description string to check.
            feat_id: Feature ID that might be in the description.

        """
        assert not fa.product_is_informative(description, feat_id)
