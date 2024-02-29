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

from contextlib import nullcontext as does_not_raise
from typing import Any, ContextManager, Dict

from deepdiff import DeepDiff
import pytest
from pytest import param

from ensembl.io.genomio.genome_metadata import dump


@pytest.mark.dependency(name="test_check_assembly_version")
@pytest.mark.parametrize(
    "genome_metadata, output, expectation",
    [
        param(
            {"assembly": {"version": "1"}},
            {"assembly": {"version": 1}},
            does_not_raise(),
            id="Version is '1'"
        ),
        param(
            {"assembly": {"accession": "GCA_00000001.1", "version": "a"}},
            {"assembly": {"accession": "GCA_00000001.1", "version": 1}},
            does_not_raise(),
            id="Version is 'a', accession is '1'",
        ),
        param(
            {"assembly": {"accession": "GCA_00000001.1"}},
            {"assembly": {"accession": "GCA_00000001.1", "version": 1}},
            does_not_raise(),
            id="No version, accession with version",
        ),
        param(
            {"assembly": {"accession": "GCA_00000001"}},
            {},
            pytest.raises(ValueError),
            id="No version, accession without version",
        ),
    ],
)
def test_check_assembly_version(
    genome_metadata: Dict[str, Any], output: Dict[str, Any], expectation: ContextManager
) -> None:
    """Tests the `dump.check_assembly_version()` method.

    Args:
        genome_metadata: Nested metadata key values from the core metadata table.
        output: Expected change in the genome metadata dictionary.
        expectation: Context manager for the expected exception (if any).
    """
    with expectation:
        dump.check_assembly_version(genome_metadata)
        assert not DeepDiff(genome_metadata, output)


@pytest.mark.dependency(name="test_check_genebuild_version")
@pytest.mark.parametrize(
    "genome_metadata, output, expectation",
    [
        param({}, {}, does_not_raise(), id="No 'genebuild' entry"),
        param(
            {"genebuild": {"version": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="Version is 'v1', no ID",
        ),
        param(
            {"genebuild": {"version": "v1", "id": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="Version is 'v1', ID dropped",
        ),
        param(
            {"genebuild": {"id": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="No version, ID moved to version",
        ),
        param({"genebuild": {}}, {}, pytest.raises(ValueError), id="No version or ID"),
    ],
)
def test_check_genebuild_version(
    genome_metadata: Dict[str, Any], output: Dict[str, Any], expectation: ContextManager
) -> None:
    """Tests the `dump.check_genebuild_version()` method.

    Args:
        genome_metadata: Nested metadata key values from the core metadata table.
        output: Expected change in the genome metadata dictionary.
        expectation: Context manager for the expected exception (if any).
    """
    with expectation:
        dump.check_genebuild_version(genome_metadata)
        assert not DeepDiff(genome_metadata, output)
