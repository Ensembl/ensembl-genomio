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
"""Unit testing of `ensembl.io.genomio.genome_metadata.dump` module.

Typical usage example::
    $ pytest test_dump.py

"""

from collections import namedtuple
from contextlib import nullcontext as does_not_raise
from typing import Any, ContextManager, Dict, List
from unittest.mock import Mock, patch

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_metadata import dump


MetaRow = namedtuple("MetaRow", "meta_key meta_value")


@pytest.mark.parametrize(
    "genome_metadata, output, expectation",
    [
        pytest.param({"assembly": {"version": "1"}}, 1, does_not_raise(), id="Version is '1'"),
        pytest.param(
            {"assembly": {"accession": "GCA_00000001.1", "version": "a"}},
            1,
            does_not_raise(),
            id="Version is 'a', accession's version is 1",
        ),
        pytest.param(
            {"assembly": {"accession": "GCA_00000001.1"}},
            1,
            does_not_raise(),
            id="No version, accession's version is 1",
        ),
        pytest.param(
            {"assembly": {"accession": "GCA_00000001"}},
            0,
            pytest.raises(ValueError),
            id="No version, accession without version",
        ),
    ],
)
def test_check_assembly_version(
    genome_metadata: Dict[str, Any], output: int, expectation: ContextManager
) -> None:
    """Tests the `dump.check_assembly_version()` method.

    Args:
        genome_metadata: Nested genome metadata key values.
        output: Expected assembly version.
        expectation: Context manager for the expected exception (if any).
    """
    with expectation:
        dump.check_assembly_version(genome_metadata)
        assert genome_metadata["assembly"]["version"] == output


@pytest.mark.parametrize(
    "genome_metadata, output, expectation",
    [
        pytest.param({}, {}, does_not_raise(), id="No 'genebuild' entry"),
        pytest.param(
            {"genebuild": {"version": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="Version is 'v1', no ID",
        ),
        pytest.param(
            {"genebuild": {"version": "v1", "id": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="Version is 'v1', ID dropped",
        ),
        pytest.param(
            {"genebuild": {"id": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="No version, ID moved to version",
        ),
        pytest.param({"genebuild": {}}, {}, pytest.raises(ValueError), id="No version or ID"),
    ],
)
def test_check_genebuild_version(
    genome_metadata: Dict[str, Any], output: Dict[str, Any], expectation: ContextManager
) -> None:
    """Tests the `dump.check_genebuild_version()` method.

    Args:
        genome_metadata: Nested genome metadata key values.
        output: Expected change in the genome metadata dictionary.
        expectation: Context manager for the expected exception (if any).
    """
    with expectation:
        dump.check_genebuild_version(genome_metadata)
        assert not DeepDiff(genome_metadata, output)


@patch("ensembl.io.genomio.genome_metadata.dump.check_genebuild_version", Mock())
@patch("ensembl.io.genomio.genome_metadata.dump.check_assembly_version", Mock())
@pytest.mark.parametrize(
    "genome_metadata, output",
    [
        ({"species": {"taxonomy_id": "5485"}}, {"species": {"taxonomy_id": 5485}}),
        ({"species": {"display_name": "Dog"}}, {"species": {"display_name": "Dog"}}),
        ({"genebuild": {"new_key": "_"}}, {"genebuild": {}}),
        ({"BRC5": "new_value"}, {}),
        ({"meta": "key", "species": {"alias": "woof"}}, {"species": {"alias": "woof"}}),
        ({"added_seq": {"region_name": [1, 2]}}, {"added_seq": {"region_name": ["1", "2"]}}),
    ],
)
def test_filter_genome_meta(genome_metadata: Dict[str, Any], output: Dict[str, Any]) -> None:
    """Tests the `dump.filter_genome_meta()` method.

    Args:
        genome_metadata: Nested genome metadata key values.
        output: Expected change in the genome metadata dictionary.
    """
    result = dump.filter_genome_meta(genome_metadata)
    assert not DeepDiff(result, output)


@patch("sqlalchemy.engine.Result")
@patch("sqlalchemy.orm.Session")
@pytest.mark.parametrize(
    "meta_data, output, expectation",
    [
        pytest.param([], {}, does_not_raise(), id="Empty meta table"),
        pytest.param(
            [
                [MetaRow("sample", "gene1")],
                [MetaRow("species.name", "dog")],
                [MetaRow("species.synonym", "puppy")],
            ],
            {"sample": "gene1", "species": {"name": "dog", "synonym": "puppy"}},
            does_not_raise(),
            id="Meta table with simple values",
        ),
        pytest.param(
            [
                [MetaRow("sample", "gene1")],
                [MetaRow("sample", "gene2")],
                [MetaRow("species.synonym", "dog")],
                [MetaRow("species.synonym", "puppy")],
            ],
            {"sample": ["gene1", "gene2"], "species": {"synonym": ["dog", "puppy"]}},
            does_not_raise(),
            id="Meta table with lists",
        ),
        pytest.param(
            [[MetaRow("species", "dog")], [MetaRow("species.synonym", "puppy")]],
            {},
            pytest.raises(ValueError),
            id="'species' and 'species.synonym' meta keys",
        ),
    ],
)
def test_get_genome_metadata(
    mock_session: Mock,
    mock_result: Mock,
    meta_data: List[MetaRow],
    output: Dict[str, Any],
    expectation: ContextManager,
) -> None:
    """Tests the `dump.get_genome_metadata()` method.

    Args:
        mock_session: A mock of `sqlalchemy.orm.Session()` class.
        meta_data: `meta` table content in a list of named tuples.
        output: Expected genome metadata dictionary.
        expectation: Context manager for the expected exception (if any).
    """
    mock_result.unique.return_value = mock_result
    mock_result.all.return_value = meta_data
    mock_session.execute.return_value = mock_result
    with expectation:
        result = dump.get_genome_metadata(mock_session)
        assert not DeepDiff(result, output)
