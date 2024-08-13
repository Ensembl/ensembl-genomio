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
"""Unit testing of `ensembl.io.genomio.seq_region.dump` module."""

from contextlib import nullcontext as no_raise
from pathlib import Path
from typing import ContextManager

import pytest
from pytest import param, raises
from sqlalchemy import text

from ensembl.io.genomio.seq_region.dump import fetch_coord_systems, fetch_seq_regions
from ensembl.core.models import CoordSystem, SeqRegion, SeqRegionAttrib, metadata
from ensembl.utils.database import UnitTestDB


@pytest.fixture(name="seq_test_db", scope="module")
def fixture_seq_test_db(db_factory) -> UnitTestDB:
    """Returns a test database with all core tables and basic data."""
    test_db: UnitTestDB = db_factory(src="no_src", name="dump_seq")
    test_db.dbc.create_all_tables(metadata)

    with test_db.dbc.session_scope() as session:
        # Add base data to dump
        if test_db.dbc.dialect == "mysql":
            session.execute(text("SET FOREIGN_KEY_CHECKS=0"))
        coord = CoordSystem(
            coord_system_id=1,
            species_id=1,
            name="primary_assembly",
            rank=1,
            attrib="default_version",
        )
        session.add(coord)
        seqr = SeqRegion(
            seq_region_id=1,
            coord_system_id=1,
            name="seqA",
            length=1000,
        )
        session.add(seqr)
    return test_db


def test_fetch_coord_systems(
    seq_test_db: UnitTestDB
) -> None:
    """Tests the `fetch_coord_system` method."""
    with seq_test_db.dbc.test_session_scope() as session:
        coords = list(fetch_coord_systems(session))
        assert len(coords) == 1
        coord: CoordSystem = coords[0]
        assert len(list(coord.seq_region)) == 1
        assert coord.name == "primary_assembly"


def test_fetch_seq_regions(
    seq_test_db: UnitTestDB
) -> None:
    """Tests the `fetch_seq_regions` method."""
    with seq_test_db.dbc.test_session_scope() as session:
        coord_system = list(fetch_coord_systems(session))[0]
        seqrs = list(fetch_seq_regions(session, coord_system))
        assert len(seqrs) == 1
        seqr = seqrs[0]
        assert seqr.name == "seqA"
