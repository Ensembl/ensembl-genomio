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
"""Unit testing of `ensembl.io.genomio.database.core_server` module."""
# pylint: disable=too-many-positional-arguments

from __future__ import annotations

from typing import Any, List, Optional

import pytest
from pytest_mock import MockerFixture
from sqlalchemy.engine import make_url

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
    """Mocker of `sqlalchemy.engine.Result` class."""

    def __init__(self, core_dbs: List[str]):
        self.core_dbs = core_dbs

    def fetchall(self) -> List[List[str]]:
        """Return a list of lists, each one containing a single core db."""
        return [[x] for x in self.core_dbs]


class MockConnection:
    """Mock a SQLAlchemy connection."""

    def __init__(self, result: MockResult) -> None:
        self.result = result

    def execute(self, *args: Any, **kwargs: Any) -> MockResult:  # pylint: disable=unused-argument
        """Returns a `MockResult` object."""
        return self.result

    def __enter__(self, *args: Any, **kwargs: Any) -> MockConnection:  # pylint: disable=unused-argument
        return self

    def __exit__(self, *args: Any, **kwargs: Any) -> None:  # pylint: disable=unused-argument
        pass


class MockEngine:
    """Mocker of `sqlalchemy.engine.Engine` class."""

    def __init__(self, core_dbs: List[str]) -> None:
        self.result = MockResult(core_dbs)

    def connect(self) -> MockConnection:
        """Return a mock connection."""
        return MockConnection(self.result)


class TestCoreServer:
    """Tests for the `CoreServer` class."""

    @pytest.mark.parametrize(
        "dbs, prefix, build, version, dbname_re, db_list, output",
        [
            ([], "", "", None, None, [], []),
            (TEST_CORES, "", None, None, "", [], TEST_CORES),
            (TEST_CORES, "prefix", "61", "111", r"_sp_", [], ["prefix_speciesF_genus_sp_core_61_111_1"]),
            (TEST_CORES, "speciesC", None, None, "", [], []),
            (TEST_CORES, "", 59, None, "", [], []),
            (TEST_CORES, "", None, 109, "", [], []),
            (TEST_CORES, "", None, None, r"_compara_", [], []),
            (TEST_CORES, "", None, None, r"", TEST_CORES[0:2], TEST_CORES[0:2]),
            (TEST_CORES, "", None, None, r"", ["nonexistent_species"], []),
        ],
    )
    def test_get_cores(
        self,
        mocker: MockerFixture,
        dbs: List[str],
        prefix: str,
        build: Optional[int],
        version: Optional[int],
        dbname_re: str,
        db_list: List[str],
        output: List[str],
    ) -> None:
        """Tests the `CoreServer.get_cores()` method.

        Args:
            mocker: Fixture to mock the connection to the server.
            dbs: Mock list of databases found in the server.
            prefix: Filter by prefix.
            build: Filter by VEuPathDB build number.
            version: Filter by Ensembl version.
            dbname_re: Filter by dbname regular expression.
            db_list: Explicit list of database names.
            output: Expected list of databases.

        """
        # Mock the engine creation that connects to the server
        mocked_engine = mocker.patch("sqlalchemy.create_engine")
        mocked_engine.return_value = MockEngine(dbs)

        # Fake server with mock get_all_core_names()
        server_url = make_url("sqlite:///:memory:")
        server = CoreServer(server_url)

        # Checks the filters from get_cores
        all_cores = server.get_cores(
            prefix=prefix, build=build, version=version, dbname_re=dbname_re, db_list=db_list
        )
        assert set(all_cores) == set(output)
