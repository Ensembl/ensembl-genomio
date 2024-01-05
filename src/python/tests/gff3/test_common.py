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
"""Unit testing of `ensembl.io.genomio.gff3.common` module."""

from typing import Optional
import pytest

from ensembl.io.genomio.gff3 import GFFMeta


class TestGFF3Merge:
    """Tests for the GFF3GeneMerger module."""

    @pytest.mark.parametrize(
        "gene_type, supported, expected_contain",
        [
            ("gene", None, "pseudogene"),
            ("gene", True, "pseudogene"),
            ("gene", False, "gap"),
            ("non_gene", True, "transposable_element"),
            ("transcript", True, "mRNA"),
            ("transcript", False, "3'UTR"),
        ],
    )
    def test_get_biotypes(self, gene_type: str, supported: Optional[bool], expected_contain: str) -> None:
        """Tests get_biotypes method."""
        if supported is None:
            result = GFFMeta.get_biotypes(gene_type)
        else:
            result = GFFMeta.get_biotypes(gene_type, supported=supported)
        assert expected_contain in result