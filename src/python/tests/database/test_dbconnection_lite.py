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
"""Unit testing of `ensembl.io.genomio.database.dbconnection_lite` module.
"""

from typing import List

from unittest.mock import patch
from sqlalchemy.engine import URL

from ensembl.io.genomio.database import DBConnectionLite


class MockURL:
    """Mocker of URL class."""

    def __init__(self, db_name: str):
        self.database = db_name


class MockEngine:
    """Mocker of `sqlalchemy.engine.Engine` class."""

    def __init__(self, db_name: str) -> None:
        self.url = MockURL(db_name)


@patch("ensembl.io.genomio.database.dbconnection_lite.create_engine")
def test_db_name(mock_create_engine) -> None:
    db_name = "testdb"
    mock_create_engine.return_value = MockEngine(db_name)

    test_url = URL.create("sqlite")
    dbc = DBConnectionLite(test_url)
    assert dbc.db_name == db_name
