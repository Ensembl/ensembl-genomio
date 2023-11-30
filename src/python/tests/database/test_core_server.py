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
"""Unit testing of :mod:`ensembl.io.genomio.database.core_server` module.

"""

from typing import List
from unittest.mock import Mock

import pytest
from pytest_mock import MockerFixture

from ensembl.io.genomio.database import CoreServer


TEST_CORES = [
    "speciesA_genus_core_60_110_1",
    "speciesB_genus_sp_core_60_110_1",
    "prefix_speciesC_genus_core_60_110_1",
    "prefix_speciesD_genus_core_61_111_1",
    "prefix_speciesE_genus_sp_core_61_110_1",
    "prefix_speciesF_genus_sp_core_61_111_1",
]


class MockResult:
    def __init__(self, core_dbs: List[str]):
        self.core_dbs = core_dbs

    def fetchall(self) -> List[List[str]]:
        return [[x] for x in self.core_dbs]


class MockEngine:
    def __init__(self, core_dbs: List[str], *args, **kwargs) -> None:
        self.result = MockResult(core_dbs)

    def execute(self, *args, **kwargs) -> MockResult:
        return self.result


class TestCoreServer:
    """Tests for the database.core_server module."""

    @pytest.mark.parametrize(
        "dbs, prefix, build, version, dbname_re, output",
        [
            ([], "", "", "", "", []),
            (TEST_CORES, "", "", "", "", TEST_CORES),
            (TEST_CORES, "prefix", "61", "111", r"_sp_", ["prefix_speciesF_genus_sp_core_61_111_1"]),
        ],
    )
    def test_get_cores(
        self,
        mocker: MockerFixture,
        dbs: List[str],
        prefix: str,
        build: str,
        version: str,
        dbname_re: str,
        output: List[str],
    ) -> None:
        """Test the CoreServer.get_cores() method.

        Args:
            ...

        """
        # Mock the engine creation that connects to the server
        mocked_engine = mocker.patch("sqlalchemy.create_engine")
        mocked_engine.return_value = MockEngine(dbs)

        # Fake server with mock get_all_core_names()
        server_url = "sqlite:///:memory:"
        server = CoreServer(server_url)

        # Checks the filters from get_cores
        all_cores = server.get_cores(prefix, build, version, dbname_re)
        assert set(all_cores) == set(output)
