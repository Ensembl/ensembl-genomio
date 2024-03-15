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

from pathlib import Path
from typing import Optional
from unittest.mock import patch

import pytest
from sqlalchemy import create_engine
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from ensembl.io.genomio.database import DBConnectionLite
from ensembl.core.models import Base, Meta

_METADATA_CONTENT = {
    "species.scientific_name": ["Lorem Ipsum"],
    "species.taxonomy_id": ["123456"],
    "species.classification": ["Insecta", "Lorem"],
}


# Create a SQLite database fixture with only a meta table and limited data
@pytest.fixture(name="db_file", scope="session")
def db_file_test(tmp_path: Path) -> Path:
    """Get a path to a db file."""
    test_db_file = tmp_path / "tmp_sqlite_core.db"
    return test_db_file


@pytest.fixture(name="db_engine", scope="session")
def db_engine_test(db_file: Path):
    """Get a SQLalchemy engine to a populated database."""
    db_url = f"sqlite:///{db_file}"
    test_db_engine = create_engine(db_url)
    Base.metadata.tables["meta"].create(test_db_engine)

    # Add some basic data
    with Session(test_db_engine) as session:
        metas = []
        for meta_key, meta_values in _METADATA_CONTENT.items():
            for meta_value in meta_values:
                metas.append(Meta(meta_key=meta_key, meta_value=meta_value))
        session.add_all(metas)
        session.commit()
    return test_db_engine


@pytest.fixture(name="dbc", scope="session")
@patch("ensembl.io.genomio.database.dbconnection_lite.create_engine")
def dbc_test(mock_create_engine, db_file: Path, db_engine) -> DBConnectionLite:
    """Provide a DBConnectionLite to the test database."""
    mock_create_engine.return_value = db_engine
    test_url = URL.create(f"sqlite://{db_file}")
    dbc = DBConnectionLite(test_url)
    return dbc


# Tests start here


def test_db_name(dbc: DBConnectionLite, db_file: Path) -> None:
    """Tests the propery db_name"""
    assert Path(dbc.db_name) == Path(db_file)


def test_get_metadata(dbc: DBConnectionLite) -> None:
    """Tests the method get_metadata()"""
    assert dbc.get_metadata() == _METADATA_CONTENT


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
def test_get_meta_value(dbc: DBConnectionLite, meta_key: str, meta_value: Optional[str]) -> None:
    """Tests the method get_meta_value()"""
    assert dbc.get_meta_value(meta_key) == meta_value


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
