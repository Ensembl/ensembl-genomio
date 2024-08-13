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
from typing import Any, ContextManager

import pytest
from pytest import param, raises
from sqlalchemy import text

from ensembl.io.genomio.seq_region.dump import (
    fetch_coord_systems,
    fetch_seq_regions,
    get_synonyms,
    get_seq_regions,
)
from ensembl.io.genomio.seq_region.external_db_map import get_external_db_map, DEFAULT_EXTERNAL_DB_MAP
from ensembl.core.models import (
    metadata,
    CoordSystem,
    SeqRegion,
    SeqRegionAttrib,
    SeqRegionSynonym,
    ExternalDb,
    AttribType,
)
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
        session.commit()

        if test_db.dbc.dialect == "mysql":
            session.execute(text("SET FOREIGN_KEY_CHECKS=1"))
    return test_db


def _add_test_synonym(session, dialect: str, synonym: str, db_name: str | None, db_id: int) -> None:
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=0"))

    dbxref_id = None
    if db_name:
        dbxref = ExternalDb(
            external_db_id=db_id,
            db_name=db_name,
            status="XREF",
            priority=1,
            db_type="PRIMARY_DB_SYNONYM",
        )
        session.add(dbxref)
        dbxref_id = 1
    seq_syn = SeqRegionSynonym(
        seq_region_id=1,
        synonym=synonym,
        external_db_id=dbxref_id,
    )
    session.add(seq_syn)
    session.commit()

    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=1"))


def _add_test_attrib(session, dialect, logic_name: str, value: str | int, attrib_id: int) -> None:
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=0"))
    attrib = AttribType(attrib_type_id=attrib_id, code=logic_name)
    session.add(attrib)
    seq_attrib = SeqRegionAttrib(attrib_type_id=attrib_id, seq_region_id=1, value=str(value))
    session.add(seq_attrib)
    session.commit()
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=1"))


@pytest.fixture(name="db_map", scope="module")
def fixture_db_map(db_factory) -> dict[str, str]:
    """Returns the default db_map."""
    return get_external_db_map(DEFAULT_EXTERNAL_DB_MAP)


def test_fetch_coord_systems(seq_test_db: UnitTestDB) -> None:
    """Tests the `fetch_coord_system` method."""
    with seq_test_db.dbc.test_session_scope() as session:
        coords = list(fetch_coord_systems(session))
        assert len(coords) == 1
        coord: CoordSystem = coords[0]
        assert len(list(coord.seq_region)) == 1
        assert coord.name == "primary_assembly"


def test_fetch_seq_regions(seq_test_db: UnitTestDB) -> None:
    """Tests the `fetch_seq_regions` method."""
    with seq_test_db.dbc.test_session_scope() as session:
        coord_system = list(fetch_coord_systems(session))[0]
        seqrs = list(fetch_seq_regions(session, coord_system))
        assert len(seqrs) == 1
        seqr = seqrs[0]
        assert seqr.name == "seqA"


def test_get_synonyms(seq_test_db: UnitTestDB, db_map: dict[str, str]) -> None:
    """Tests the `get_synonyms` method."""
    with seq_test_db.dbc.test_session_scope() as session:
        _add_test_synonym(session, seq_test_db.dbc.dialect, "synA", "db1", 1)
        coord_system = list(fetch_coord_systems(session))[0]
        seqr = list(fetch_seq_regions(session, coord_system))[0]
        syns = get_synonyms(seqr, db_map)
        assert len(syns) == 1
        syn = syns[0]
        assert syn == {"name": "synA", "source": "db1"}


@pytest.mark.parametrize(
    "attribs, expected_output",
    [
        param({}, None, id="no toplevel"),
        param({"codon_table": 1}, None, id="no toplevel, other attribs"),
        param({"toplevel": 1}, {}, id="1 toplevel"),
        param({"toplevel": 1, "codon_table": 2}, {"codon_table": 2}, id="With codon table"),
        param({"toplevel": 1, "sequence_location": "chr"}, {"location": "chr"}, id="Location"),
        param({"toplevel": 1, "circular_seq": 1}, {"circular": 1}, id="circular"),
        param({"toplevel": 1, "circular_seq": 0}, {}, id="not circular"),
        param({"toplevel": 1, "coord_system_tag": "contig"}, {"coord_system_level": "contig"}, id="contig level"),
    ],
)
def test_get_seq_regions_attribs(
    seq_test_db: UnitTestDB,
    db_map: dict[str, str],
    attribs: dict[str, str],
    expected_output: dict[str, Any] | None,
) -> None:
    """Tests the `get_seq_regions` method."""
    expected_dict = {"name": "seqA", "length": 1000, "coord_system_level": "primary_assembly"}
    with seq_test_db.dbc.test_session_scope() as session:
        attrib_num = 1
        for attrib_name, attrib_value in attribs.items():
            _add_test_attrib(session, seq_test_db.dbc.dialect, attrib_name, attrib_value, attrib_num)
            attrib_num += 1
        output_seqr = get_seq_regions(session, db_map)
        if expected_output is None:
            assert output_seqr == []
        else:
            expected_dict.update(expected_output)
            assert output_seqr[0] == expected_dict


def test_get_seq_regions(
    seq_test_db: UnitTestDB,
    db_map: dict[str, str],
) -> None:
    """Tests the `get_seq_regions` method."""
    expected_dict = {
        "name": "seqA",
        "length": 1000,
        "coord_system_level": "primary_assembly",
        "synonyms": [{"name": "synA", "source": "GenBank"}, {"name": "synB"}],
        "added_sequence": {"accession": "ADDED1"},
    }
    with seq_test_db.dbc.test_session_scope() as session:
        _add_test_attrib(session, seq_test_db.dbc.dialect, "toplevel", 1, 1)
        _add_test_attrib(session, seq_test_db.dbc.dialect, "added_seq_accession", "ADDED1", 2)
        _add_test_synonym(session, seq_test_db.dbc.dialect, "synA", "GenBank", 1)
        _add_test_synonym(session, seq_test_db.dbc.dialect, "synB", None, 2)
        output_seqr = get_seq_regions(session, db_map)
        assert output_seqr[0] == expected_dict
