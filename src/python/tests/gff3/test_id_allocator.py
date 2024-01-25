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
"""Unit testing of `ensembl.io.genomio.gff3.id_allocator` module."""

from typing import List

import pytest

from ensembl.io.genomio.gff3.id_allocator import IDAllocator


@pytest.mark.parametrize(
    "prefix, expected_ids",
    [
        pytest.param(None, ["TMP_1", "TMP_2"], id="Default prefix"),
        pytest.param("MYPREF_", ["MYPREF_1", "MYPREF_2"], id="Prefix MYPREF_"),
    ],
)
def test_generate_id(prefix: str, expected_ids: List[str]) -> None:
    """Test IDs generation."""
    ids = IDAllocator()
    if prefix:
        ids.prefix = prefix

    id1 = ids.generate_id()
    id2 = ids.generate_id()
    new_ids = [id1, id2]
    assert new_ids == expected_ids

@pytest.mark.parametrize(
    "min_id_length, test_id, outcome",
    [
        pytest.param(None, "LOREMIPSUM_01", True, id="OK ID"),
        pytest.param(None, "", False, id="Empty"),
        pytest.param(None, "A", False, id="Too short"),
        pytest.param(None, "Abc", False, id="Too short 2"),
        pytest.param(5, "Abcde", True, id="At custom min length"),
        pytest.param(5, "Abcd", False, id="Below custom min length"),
        pytest.param(None, "CHR1:100..200", False, id="Coordinates"),
        pytest.param(None, "LOREM|IPSUM", False, id="Special char |"),
        pytest.param(None, "LOREM IPSUM", False, id="Special char space"),
        pytest.param(None, "Trnaa-UAA", False, id="Trna ID, lower case"),
        pytest.param(None, "TRNAA-UAA", False, id="Trna ID, upper case"),
    ],
)
def test_valid_id(min_id_length: int, test_id: str, outcome: bool) -> None:
    """Test ID validity check."""
    ids = IDAllocator()
    if min_id_length:
        ids.min_id_length = min_id_length

    assert ids.valid_id(test_id) == outcome

@pytest.mark.parametrize(
    "test_id, prefixes, outcome",
    [
        pytest.param("LOREM-IPSUM1", [], "LOREM-IPSUM1", id="No prefixes"),
        pytest.param("LOREM-IPSUM1", ["DOLOR"], "LOREM-IPSUM1", id="Unused prefix"),
        pytest.param("LOREM-IPSUM1", ["LOREM-"], "IPSUM1", id="Found 1 prefix"),
        pytest.param("LOREM-IPSUM1", ["LOREM-", "IPSUM"], "IPSUM1", id="Only 1 prefix is removed"),
    ],
)
def test_remove_prefixes(test_id: int, prefixes: List[str], outcome: str) -> None:
    """Test prefix removal."""
    ids = IDAllocator()

    assert ids.remove_prefixes(test_id, prefixes) == outcome
