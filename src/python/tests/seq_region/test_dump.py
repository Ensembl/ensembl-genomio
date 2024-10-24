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
# pylint: disable=too-many-positional-arguments

from typing import Any, Callable

import pytest
from pytest import param
from sqlalchemy import text
from sqlalchemy.orm import Session

from ensembl.core.models import (
    metadata,
    AttribType,
    CoordSystem,
    ExternalDb,
    Karyotype,
    SeqRegion,
    SeqRegionAttrib,
    SeqRegionSynonym,
)
from ensembl.io.genomio.seq_region.dump import (
    fetch_coord_systems,
    fetch_seq_regions,
    get_added_sequence,
    get_karyotype,
    get_seq_regions,
    get_synonyms,
)
from ensembl.io.genomio.external_db.db_map import DEFAULT_EXTERNAL_DB_MAP, get_external_db_map
from ensembl.utils.database import UnitTestDB


@pytest.fixture(name="seq_test_db", scope="module")
def fixture_seq_test_db(db_factory: Callable) -> UnitTestDB:
    """Returns a test database with all core tables and basic data (1 coord_system + 1 seq_region)."""
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


def _add_test_synonym(session: Session, dialect: str, synonym: str, db_name: str | None, db_id: int) -> None:
    """Add a seq_region synonym to the test seq_region.

    Args:
        session: SQLalchemy session.
        dialect: Dialect of the DB connection.
        synonym: Synonym name.
        db_name: DB xref name.
        db_id: DB xref ID (must be unique).

    """
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


def _add_test_attrib(
    session: Session, dialect: str, logic_name: str, value: str | int, attrib_id: int
) -> None:
    """Add a seq_region attrib to the test seq_region.

    Args:
        session: SQLalchemy session.
        dialect: Dialect of the DB connection.
        logic_name: Attrib type logic name.
        value: Attrib value.
        attrib_id: Attrib ID (must be unique).

    """
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=0"))
    attrib = AttribType(attrib_type_id=attrib_id, code=logic_name)
    session.add(attrib)
    seq_attrib = SeqRegionAttrib(attrib_type_id=attrib_id, seq_region_id=1, value=str(value))
    session.add(seq_attrib)
    session.commit()
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=1"))


def _add_test_karyotype(
    session: Session, dialect: str, start: int, end: int, band: str | None = None, stain: str | None = None
) -> None:
    """Add a seq_region karyotype band to the test seq_region.

    Args:
        session: SQLalchemy session.
        dialect: Dialect of the DB connection.
        start: Start of the band.
        end: End of the band.
        band: Name of the band.
        stain: Name of the stain.

    """
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=0"))
    karyo = Karyotype(seq_region_id=1, seq_region_start=start, seq_region_end=end, band=band, stain=stain)
    session.add(karyo)
    session.commit()
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=1"))


@pytest.fixture(name="db_map", scope="module")
def fixture_db_map() -> dict[str, str]:
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
    "bands, expected_kar",
    [
        param([], [], id="no karyotype"),
        param([{}], [{"start": 1, "end": 10}], id="no name/stain"),
        param([{"name": "bandA"}], [{"start": 1, "end": 10, "name": "bandA"}], id="one band"),
        param([{"stain": "stainA"}], [{"start": 1, "end": 10, "stain": "stainA"}], id="one stain"),
        param(
            [{"stain": "TEL"}],
            [{"start": 1, "end": 10, "stain": "TEL", "structure": "telomere"}],
            id="telomere",
        ),
    ],
)
def test_get_karyotype(seq_test_db: UnitTestDB, bands: list, expected_kar: dict) -> None:
    """Tests the `get_synonyms` method.

    Args:
        bands: List of dicts with keys `name`, `stain` for each test band to add.
        expected_kar: List of expected karyotype dicts output.

    """
    with seq_test_db.dbc.test_session_scope() as session:
        for band in bands:
            _add_test_karyotype(session, seq_test_db.dbc.dialect, 1, 10, band.get("name"), band.get("stain"))
        coord_system = list(fetch_coord_systems(session))[0]
        seqr = list(fetch_seq_regions(session, coord_system))[0]
        kar = get_karyotype(seqr)
        assert kar == expected_kar


@pytest.mark.parametrize(
    "attribs, expected_output",
    [
        param({}, None, id="no toplevel"),
        param({"codon_table": 1}, None, id="no toplevel, other attribs"),
        param({"toplevel": 1}, {}, id="1 toplevel"),
        param({"toplevel": 1, "codon_table": 2}, {"codon_table": 2}, id="with codon table"),
        param({"toplevel": 1, "sequence_location": "chr"}, {"location": "chr"}, id="location"),
        param({"toplevel": 1, "circular_seq": 1}, {"circular": 1}, id="circular"),
        param({"toplevel": 1, "circular_seq": 0}, {}, id="not circular"),
        param(
            {"toplevel": 1, "coord_system_tag": "contig"}, {"coord_system_level": "contig"}, id="contig level"
        ),
    ],
)
def test_get_seq_regions_attribs(
    seq_test_db: UnitTestDB,
    db_map: dict[str, str],
    attribs: dict[str, str],
    expected_output: dict[str, Any] | None,
) -> None:
    """Tests the `get_seq_regions_attribs` method.

    Args:
        attribs: Dicts with all attrib-value pairs for each attrib to add to the seq_region.
        expected_kar: List of expected karyotype dicts output.

    """
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


@pytest.mark.parametrize(
    "attribs, expected_output",
    [
        param({}, {}, id="no added seq"),
        param({"added_seq_accession": "ACC1"}, {"accession": "ACC1"}, id="added accession"),
        param(
            {
                "added_seq_accession": "ACC1",
                "added_seq_asm_pr_nam": "asm_name",
                "added_seq_asm_pr_url": "asm_url",
            },
            {"accession": "ACC1", "assembly_provider": {"name": "asm_name", "url": "asm_url"}},
            id="added accession + assembly provider",
        ),
        param(
            {
                "added_seq_accession": "ACC1",
                "added_seq_ann_pr_nam": "ann_name",
                "added_seq_ann_pr_url": "ann_url",
            },
            {"accession": "ACC1", "annotation_provider": {"name": "ann_name", "url": "ann_url"}},
            id="added accession + annotation provider",
        ),
        param(
            {
                "added_seq_accession": "ACC1",
                "added_seq_asm_pr_nam": "asm_name",
                "added_seq_ann_pr_url": "ann_url",
            },
            {"accession": "ACC1"},
            id="added accession + incomplete annotation/assembly provider",
        ),
    ],
)
def test_get_added_sequence(
    seq_test_db: UnitTestDB, attribs: dict[str, str], expected_output: dict[str, str | dict]
) -> None:
    """Tests the `get_added_sequences` method.

    Args:
        attribs: Dicts with all attrib-value pairs for each attrib to add to the seq_region.
        expected_output: Expected output.

    """
    with seq_test_db.dbc.test_session_scope() as session:
        attrib_num = 1
        for attrib_name, attrib_value in attribs.items():
            _add_test_attrib(session, seq_test_db.dbc.dialect, attrib_name, attrib_value, attrib_num)
            attrib_num += 1
        coord_system = list(fetch_coord_systems(session))[0]
        seqr = list(fetch_seq_regions(session, coord_system))[0]
        added = get_added_sequence(seqr)
        assert added == expected_output


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
        "karyotype_bands": [{"start": 1, "end": 10, "name": "bandA", "stain": "stainA"}],
    }
    with seq_test_db.dbc.test_session_scope() as session:
        _add_test_attrib(session, seq_test_db.dbc.dialect, "toplevel", 1, 1)
        _add_test_attrib(session, seq_test_db.dbc.dialect, "added_seq_accession", "ADDED1", 2)
        _add_test_synonym(session, seq_test_db.dbc.dialect, "synA", "GenBank", 1)
        _add_test_synonym(session, seq_test_db.dbc.dialect, "synB", None, 2)
        _add_test_karyotype(session, seq_test_db.dbc.dialect, 1, 10, "bandA", "stainA")
        output_seqr = get_seq_regions(session, db_map)
        assert output_seqr[0] == expected_dict
