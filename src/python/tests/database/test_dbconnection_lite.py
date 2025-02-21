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
"""Unit testing of `ensembl.io.genomio.database.dbconnection_lite` module."""

from typing import Callable, Optional

import pytest

from ensembl.io.genomio.database import DBConnectionLite
from ensembl.core.models import CoordSystem, Meta, metadata
from ensembl.utils.database import UnitTestDB


_METADATA_CONTENT = {
    "species.scientific_name": ["Lorem Ipsum"],
    "species.taxonomy_id": ["123456"],
    "species.classification": ["Insecta", "Lorem"],
}


@pytest.fixture(name="meta_test_db", scope="module")
def fixture_meta_test_db(db_factory: Callable) -> UnitTestDB:
    """Returns a test database with a meta table and basic data."""
    test_db: UnitTestDB = db_factory("", "get_metadata")
    test_db.dbc.create_table(metadata.tables["coord_system"])
    test_db.dbc.create_table(metadata.tables["meta"])
    # Add data
    with test_db.dbc.session_scope() as session:
        session.add(CoordSystem(species_id=1, name="Foo", rank=1))
        for meta_key, meta_values in _METADATA_CONTENT.items():
            for meta_value in meta_values:
                session.add(Meta(species_id=1, meta_key=meta_key, meta_value=meta_value))
        session.commit()
    return test_db


# Use ensembl-utils UnitTestDB
def test_get_metadata(meta_test_db: UnitTestDB) -> None:
    """Tests the method get_metadata()"""

    # Check the new connection lite
    dblite = DBConnectionLite(meta_test_db.dbc.url)
    assert dblite.get_metadata() == _METADATA_CONTENT


@pytest.mark.parametrize(
    "meta_key, meta_value",
    [
        pytest.param(
            "species.scientific_name", _METADATA_CONTENT["species.scientific_name"][0], id="Unique key exists"
        ),
        pytest.param(
            "species.classification", _METADATA_CONTENT["species.classification"][0], id="First key exists"
        ),
        pytest.param("lorem.ipsum", None, id="Non-existing key, 2 parts"),
        pytest.param("lorem_ipsum", None, id="Non-existing key, 1 part"),
    ],
)
def test_get_meta_value(meta_test_db: UnitTestDB, meta_key: str, meta_value: Optional[str]) -> None:
    """Tests the method get_meta_value()"""
    dblite = DBConnectionLite(meta_test_db.dbc.url)
    assert dblite.get_meta_value(meta_key) == meta_value


@pytest.mark.parametrize(
    "db_name, release_version",
    [
        pytest.param("coredb_core_66_111_1", "66", id="Release version in name"),
        pytest.param("coredb_core_111_1", "", id="No release version in name"),
        pytest.param("lorem_ipsum_66_111_1", "", id="Some release version but wrong format for the rest"),
    ],
)
def test_get_project_release(db_name: str, release_version: str) -> None:
    """Tests the method get_project_release()."""
    db_url = f"sqlite:///{db_name}"
    dbc = DBConnectionLite(db_url)
    assert dbc.get_project_release() == release_version
