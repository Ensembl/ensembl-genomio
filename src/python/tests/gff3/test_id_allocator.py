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
