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
import pytest
from unittest.mock import patch

from sqlalchemy import create_engine
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from ensembl.io.genomio.database import DBConnectionLite
from ensembl.core.models import Base, Meta

metadata_content = {
    "species": ["Lorem Ipsum"],
    "classification": ["Insecta", "Lorem"],
}

# Create a database with only a meta table
@pytest.fixture(name="db_file", scope='session')
def db_file_test(tmp_dir: Path) -> Path:
    """Get a path to a db file."""
    test_db_file = tmp_dir / "tmp_sqlite.db"
    return test_db_file

@pytest.fixture(name="db_engine", scope='session')
def db_engine_test(db_file: Path):
    """Get a SQLalchemy engine to a populated database."""
    db_url = f"sqlite:///{db_file}"
    test_db_engine = create_engine(db_url)
    Base.metadata.tables["meta"].create(test_db_engine)

    # Add some basic data
    with Session(test_db_engine) as session:
        metas = []
        for meta_key, meta_values in metadata_content.items():
            for meta_value in meta_values:
                metas.append(Meta(meta_key=meta_key, meta_value=meta_value))
        session.add_all(metas)
        session.commit()

    return test_db_engine


@pytest.fixture(name="dbc", scope="session")
@patch("ensembl.io.genomio.database.dbconnection_lite.create_engine")
def dbc_test(mock_create_engine, db_file: Path, db_engine) -> None:
    mock_create_engine.return_value = db_engine
    test_url = URL.create(f"sqlite://{db_file}")
    dbc = DBConnectionLite(test_url)
    return dbc

def test_db_name(dbc, db_file: Path) -> None:
    assert Path(dbc.db_name) == Path(db_file)


def test_load_metadata(dbc) -> None:
    dbc._load_metadata()
    assert dbc._metadata == metadata_content
