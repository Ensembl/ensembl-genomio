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
"""Unit testing of the GenomIO feature database loader."""

# This module keeps mapping, persistence, and orchestration coverage together
# because they exercise one loader contract; the size warning is intentional.
# pylint: disable=too-many-lines

from collections.abc import Generator, Mapping
from contextlib import contextmanager
from datetime import datetime
import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace
from typing import cast

import pytest
from sqlalchemy import create_engine, make_url, select
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from ensembl.core.models import (
    Analysis,
    AnalysisDescription,
    AttribType,
    DNAAlignFeatureAttrib,
    DnaAlignFeature,
    RepeatConsensus,
    RepeatFeature,
    SeqRegion,
)
from ensembl.io.genomio.features.load import core_load

JsonValue = core_load.JsonValue


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param("+", 1, id="positive strand"),
        pytest.param("-", -1, id="negative strand"),
        pytest.param(".", 0, id="no strand")
    ],
)
def test_normalise_strand(value: str, expected: int) -> None:
    """GenomIO strand values map to Ensembl integer values.

    Args:
        value: GenomIO strand symbol.
        expected: Expected Ensembl integer strand.
    """
    assert core_load._normalise_strand(value) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param("text", "text", id="string value"),
        pytest.param(True, "true", id="boolean true"),
        pytest.param(False, "false", id="boolean false"),
        pytest.param(1, "1", id="integer value"),
        pytest.param(2.5, "2.5", id="float value"),
    ],
)
def test_attribute_scalar_to_text(value: object, expected: str) -> None:
    """Scalar ncRNA attributes are serialized as database text.

    Args:
        value: Input scalar.
        expected: Expected serialized text.
    """
    assert core_load._attribute_scalar_to_text(value) == expected


def test_ncrna_row_uses_cmscan_model_coordinates_and_attributes() -> None:
    """cmscan fields map to dna_align_feature and its attribute table."""
    feature: dict[str, JsonValue] = {
        "seq_region": "chr1",
        "seq_region_start": 10,
        "seq_region_end": 30,
        "seq_region_strand": "-",
        "biotype": "miRNA",
        "display_label": "Rfam test",
        "score": 86.9,
        "evalue": 2.2e-15,
        "target_name": "RF00001",
        "model_start": 3,
        "model_end": 20,
        "is_significant": True,
    }

    row, attributes = core_load._ncrna_row(feature, "cmscan", 9, 7)

    assert row == {
        "seq_region_id": 9,
        "seq_region_start": 10,
        "seq_region_end": 30,
        "seq_region_strand": -1,
        "hit_start": 3,
        "hit_end": 20,
        "hit_strand": 1,
        "hit_name": "RF00001",
        "analysis_id": 7,
        "score": 86.9,
        "evalue": 2.2e-15,
        "perc_ident": None,
        "cigar_line": None,
        "external_db_id": None,
        "hcoverage": None,
        "align_type": "ensembl",
    }
    assert attributes == {
        "biotype": "miRNA",
        "display_label": "Rfam test",
        "is_significant": "true",
        "ncrna_tool": "cmscan",
    }


@pytest.mark.parametrize(
    ("tool", "feature", "expected_name"),
    [
        pytest.param(
            "trnascan",
            {"display_label": "tRNA-Phe", "isotype": "Phe"},
            "tRNA-Phe",
            id="trnascan with display label"
        ),
        pytest.param(
            "trnascan",
            {"isotype": "Phe"},
            "Phe",
            id="trnascan with isotype"
        ),
        pytest.param(
            "mirbase",
            {"mirna_name": "ath-miR156a", "mirna_type": "hairpin"},
            "ath-miR156a",
            id="mirbase with mirna name"
        ),
    ],
)
def test_ncrna_row_uses_tool_specific_name_and_feature_length(
    tool: str, feature: dict[str, str], expected_name: str
) -> None:
    """Non-cmscan ncRNA records use feature length as hit coordinates.

    Args:
        tool: ncRNA tool under test.
        feature: Tool-specific feature fields.
        expected_name: Expected alignment hit name.
    """
    base: dict[str, JsonValue] = {
        "seq_region": "chr1",
        "seq_region_start": 10,
        "seq_region_end": 30,
        "seq_region_strand": "+",
        "biotype": "ncRNA",
        **feature,
    }
    if tool == "trnascan":
        base.update({"anticodon": "GAA"})
    else:
        base.update({"mirna_type": "hairpin"})

    row, attributes = core_load._ncrna_row(base, tool, 9, 7)

    assert row["hit_name"] == expected_name
    assert row["hit_start"] == 1
    assert row["hit_end"] == 21
    assert attributes["ncrna_tool"] == tool


def test_ncrna_row_rejects_hit_name_over_40_characters() -> None:
    """The core hit_name column length is enforced without truncation."""
    feature: dict[str, JsonValue] = {
        "seq_region": "chr1",
        "seq_region_start": 1,
        "seq_region_end": 10,
        "seq_region_strand": "+",
        "biotype": "ncRNA",
        "target_name": "x" * 41,
        "is_significant": True,
    }

    with pytest.raises(ValueError, match="at most 40"):
        core_load._ncrna_row(feature, "cmscan", 9, 7)


def test_parse_run_date_normalises_z_and_offset_values() -> None:
    """Run dates are converted to UTC-naive datetimes."""
    assert core_load._parse_run_date("2026-02-18T00:00:00Z") == datetime(2026, 2, 18)
    assert core_load._parse_run_date("2026-02-18T01:00:00+01:00") == datetime(2026, 2, 18)
    assert core_load._parse_run_date("2026-02-18T00:00:00") == datetime(2026, 2, 18)


def test_synthetic_consensus_uses_analysis_logic_name() -> None:
    """Missing repeat consensus gets the agreed synthetic values."""
    assert core_load._synthetic_consensus("repeat_logic") == {
        "repeat_name": "repeat_logic",
        "repeat_class": "repeat_logic",
        "repeat_type": "repeat_logic",
        "repeat_consensus": "N",
    }


def test_repeat_row_maps_core_columns() -> None:
    """Repeat fields map to the repeat_feature columns."""
    feature: dict[str, JsonValue] = {
        "seq_region": "chr1",
        "seq_region_start": 1,
        "seq_region_end": 10,
        "seq_region_strand": ".",
        "repeat_start": 2,
        "repeat_end": 9,
        "score": 4.5,
    }

    assert core_load._repeat_row(feature, 7, 9, 11) == {
        "seq_region_id": 9,
        "seq_region_start": 1,
        "seq_region_end": 10,
        "seq_region_strand": 0,
        "repeat_start": 2,
        "repeat_end": 9,
        "repeat_consensus_id": 11,
        "analysis_id": 7,
        "score": 4.5,
    }


def _write_input_file(tmp_path: Path) -> Path:
    """Create a readable input path for parser tests.

    Args:
        tmp_path: Pytest temporary directory.

    Returns:
        Path to the generated input file.
    """
    input_path = tmp_path / "features.json"
    input_path.write_text("{}", encoding="utf-8")
    return input_path


def test_parse_args_uses_required_command_line_connection_values(tmp_path: Path) -> None:
    """Parser builds the database URL from explicit command-line values.

    Args:
        tmp_path: Pytest temporary directory.
    """
    args = core_load.parse_args(
        [
            "--input_file",
            str(_write_input_file(tmp_path)),
            "--db_host",
            "cli-host",
            "--db_port",
            "3307",
            "--db_user",
            "cli-user",
            "--db_password",
            "cli-secret",
            "--db_database",
            "cli-db",
            "--batch-size",
            "25",
            "--skip-validation",
        ]
    )

    assert args.db_host == "cli-host"
    assert args.db_port == 3307
    assert args.db_user == "cli-user"
    assert args.db_password == "cli-secret"
    assert args.db_database == "cli-db"
    assert args.db_url == URL.create("mysql", "cli-user", "cli-secret", "cli-host", 3307, "cli-db")
    assert args.batch_size == 25
    assert args.skip_validation is True


@pytest.mark.parametrize(
    ("missing"),
    [
        pytest.param("--db_host", id="missing db host"),
        pytest.param("--db_port", id="missing db port"),
        pytest.param("--db_user", id="missing db user"),
        pytest.param("--db_database", id="missing db database")
    ]
)
def test_parse_args_requires_each_connection_parameter(missing: str, tmp_path: Path) -> None:
    """Each required connection option is rejected when omitted.

    Args:
        missing: Required option to omit.
        tmp_path: Pytest temporary directory.

    Raises:
        SystemExit: Expected when a required option is missing.
    """
    arguments = [
        "--input_file",
        str(_write_input_file(tmp_path)),
        "--db_host",
        "db.example",
        "--db_port",
        "3306",
        "--db_user",
        "loader",
        "--db_database",
        "species_core",
    ]
    position = arguments.index(missing)
    del arguments[position : position + 2]

    with pytest.raises(SystemExit):
        core_load.parse_args(arguments)


def test_parse_args_rejects_nonpositive_batch_size(tmp_path: Path) -> None:
    """Batch size must be positive.

    Args:
        tmp_path: Pytest temporary directory.

    Raises:
        SystemExit: Expected for a zero batch size.
    """
    with pytest.raises(SystemExit):
        core_load.parse_args(
            [
                "--input_file",
                str(_write_input_file(tmp_path)),
                "--db_host",
                "db.example",
                "--db_port",
                "3306",
                "--db_user",
                "loader",
                "--db_database",
                "species_core",
                "--batch-size",
                "0",
            ]
        )


@pytest.fixture(name="parent_session")
def fixture_parent_session() -> Session:
    """Return a SQLite session containing parent loader tables.

    Returns:
        Session connected to an isolated SQLite test database.
    """
    engine = create_engine("sqlite:///:memory:")
    with engine.begin() as connection:
        connection.exec_driver_sql(
            "CREATE TABLE analysis ("
            "analysis_id INTEGER PRIMARY KEY, created DATETIME, logic_name VARCHAR(128) NOT NULL UNIQUE, "
            "db VARCHAR(120), db_version VARCHAR(40), db_file VARCHAR(120), program VARCHAR(80), "
            "program_version VARCHAR(40), program_file VARCHAR(80), parameters TEXT, module VARCHAR(80), "
            "module_version VARCHAR(40), gff_source VARCHAR(40), gff_feature VARCHAR(40))"
        )
        connection.exec_driver_sql(
            "CREATE TABLE analysis_description ("
            "analysis_id INTEGER PRIMARY KEY, description TEXT, display_label VARCHAR(255) NOT NULL, "
            "displayable BOOLEAN DEFAULT 1, web_data TEXT, "
            "FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id))"
        )
        connection.exec_driver_sql(
            "CREATE TABLE repeat_consensus ("
            "repeat_consensus_id INTEGER PRIMARY KEY, repeat_name VARCHAR(255) NOT NULL, "
            "repeat_class VARCHAR(100) NOT NULL, repeat_type VARCHAR(40) NOT NULL, repeat_consensus TEXT)"
        )
    return Session(engine)


def _analysis_payload(logic_name: str = "repeat_logic") -> dict[str, JsonValue]:
    """Return valid analysis data for persistence tests.

    Args:
        logic_name: Analysis logic name.

    Returns:
        Valid analysis metadata.
    """
    return {
        "run_date": "2026-02-18T00:00:00Z",
        "logic_name": logic_name,
        "display_label": "Repeat test",
        "description": "Repeat description",
        "program": "RepeatMasker",
        "program_version": "4.1.7",
        "program_parameters": "-lib custom",
        "db": "library",
        "db_version": "1",
    }


def _consensus_payload(key: str = "a" * 64) -> dict[str, JsonValue]:
    """Return a consensus object for persistence tests.

    Args:
        key: Input consensus key.

    Returns:
        Valid repeat consensus metadata.
    """
    return {
        "repeat_consensus_key": key,
        "repeat_name": "Alu",
        "repeat_class": "SINE",
        "repeat_type": "Alu",
        "repeat_consensus": "ACGT",
    }


def test_get_or_create_analysis_creates_parent_rows(parent_session: Session) -> None:
    """New analysis metadata creates analysis and description rows.

    Args:
        parent_session: SQLite session containing parent tables.
    """
    analysis_id = core_load._get_or_create_analysis(parent_session, _analysis_payload())

    analysis = parent_session.get(Analysis, analysis_id)
    description = parent_session.scalar(
        select(AnalysisDescription).where(AnalysisDescription.analysis_id == analysis_id)
    )
    assert analysis is not None
    assert analysis.logic_name == "repeat_logic"
    assert analysis.parameters == "-lib custom"
    assert description is not None
    assert description.description == "Repeat description"


def test_get_or_create_analysis_reuses_matching_metadata(parent_session: Session) -> None:
    """Matching analysis metadata reuses the existing analysis ID.

    Args:
        parent_session: SQLite session containing parent tables.
    """
    first_id = core_load._get_or_create_analysis(parent_session, _analysis_payload())
    second_id = core_load._get_or_create_analysis(parent_session, _analysis_payload())

    assert second_id == first_id
    assert parent_session.query(Analysis).count() == 1
    assert parent_session.query(AnalysisDescription).count() == 1


@pytest.mark.parametrize(
    ("field"),
    [
        pytest.param("program", id="program"),
        pytest.param("description", id="description"),
        pytest.param("display_label", id="display_label")
    ]
)
def test_get_or_create_analysis_rejects_conflicting_metadata(parent_session: Session, field: str) -> None:
    """Conflicting analysis metadata is rejected rather than updated.

    Args:
        parent_session: SQLite session containing parent tables.
        field: Analysis field to alter.

    Raises:
        ValueError: Expected for conflicting metadata.
    """
    original = _analysis_payload()
    core_load._get_or_create_analysis(parent_session, original)
    conflicting = dict(original)
    conflicting[field] = f"different {field}"

    with pytest.raises(ValueError, match="conflicts"):
        core_load._get_or_create_analysis(parent_session, conflicting)


def test_get_or_create_analysis_rejects_missing_description(parent_session: Session) -> None:
    """An existing analysis without its description row is invalid.

    Args:
        parent_session: SQLite session containing parent tables.

    Raises:
        ValueError: Expected when the description row is absent.
    """
    analysis_id = core_load._get_or_create_analysis(parent_session, _analysis_payload())
    parent_session.query(AnalysisDescription).filter_by(analysis_id=analysis_id).delete()
    parent_session.flush()

    with pytest.raises(ValueError, match="no analysis_description"):
        core_load._get_or_create_analysis(parent_session, _analysis_payload())


def test_consensus_index_context_cleans_up_on_exception() -> None:
    """The temporary consensus index closes on exceptional exit."""
    connection: sqlite3.Connection | None = None

    with pytest.raises(RuntimeError, match="stop"):
        with core_load._consensus_index() as index:
            connection = index
            index.execute("SELECT 1")
            raise RuntimeError("stop")

    assert connection is not None
    with pytest.raises(sqlite3.ProgrammingError, match="closed"):
        connection.execute("SELECT 1")


def test_get_or_create_consensus_creates_and_reuses_row(parent_session: Session) -> None:
    """A consensus key maps to one reusable core row.

    Args:
        parent_session: SQLite session containing parent tables.
    """
    with core_load._consensus_index() as index:
        first_id = core_load._get_or_create_consensus(parent_session, _consensus_payload(), index)
        second_id = core_load._get_or_create_consensus(parent_session, _consensus_payload(), index)

    assert second_id == first_id
    assert parent_session.query(RepeatConsensus).count() == 1


def test_get_or_create_consensus_rejects_key_conflict(parent_session: Session) -> None:
    """A repeated key with different consensus data fails.

    Args:
        parent_session: SQLite session containing parent tables.

    Raises:
        ValueError: Expected when a key's consensus data changes.
    """
    first = _consensus_payload()
    conflicting = dict(first)
    conflicting["repeat_consensus"] = "TGCA"
    with core_load._consensus_index() as index:
        core_load._get_or_create_consensus(parent_session, first, index)

        with pytest.raises(ValueError, match="consensus key"):
            core_load._get_or_create_consensus(parent_session, conflicting, index)


def test_get_or_create_consensus_reuses_matching_database_row(parent_session: Session) -> None:
    """A matching consensus row is reused even for a new input key.

    Args:
        parent_session: SQLite session containing parent tables.
    """
    existing = RepeatConsensus(
        repeat_name="Alu",
        repeat_class="SINE",
        repeat_type="Alu",
        repeat_consensus="ACGT",
    )
    parent_session.add(existing)
    parent_session.flush()

    with core_load._consensus_index() as index:
        consensus_id = core_load._get_or_create_consensus(parent_session, _consensus_payload("b" * 64), index)

    assert consensus_id == existing.repeat_consensus_id


@pytest.fixture(name="feature_session")
def fixture_feature_session() -> Session:
    """Return a SQLite session containing the loader's feature tables.

    Returns:
        Session connected to an isolated SQLite test database.
    """
    engine = create_engine("sqlite:///:memory:")
    with engine.begin() as connection:
        connection.exec_driver_sql(
            "CREATE TABLE analysis ("
            "analysis_id INTEGER PRIMARY KEY, created DATETIME, logic_name VARCHAR(128) NOT NULL UNIQUE, "
            "db VARCHAR(120), db_version VARCHAR(40), db_file VARCHAR(120), program VARCHAR(80), "
            "program_version VARCHAR(40), program_file VARCHAR(80), parameters TEXT, module VARCHAR(80), "
            "module_version VARCHAR(40), gff_source VARCHAR(40), gff_feature VARCHAR(40))"
        )
        connection.exec_driver_sql(
            "CREATE TABLE analysis_description ("
            "analysis_id INTEGER PRIMARY KEY, description TEXT, display_label VARCHAR(255) NOT NULL, "
            "displayable BOOLEAN DEFAULT 1, web_data TEXT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE repeat_consensus ("
            "repeat_consensus_id INTEGER PRIMARY KEY, repeat_name VARCHAR(255) NOT NULL, "
            "repeat_class VARCHAR(100) NOT NULL, repeat_type VARCHAR(40) NOT NULL, repeat_consensus TEXT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE seq_region ("
            "seq_region_id INTEGER PRIMARY KEY, name VARCHAR(255) NOT NULL, "
            "coord_system_id INTEGER NOT NULL, length INTEGER NOT NULL)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE repeat_feature ("
            "repeat_feature_id INTEGER PRIMARY KEY, seq_region_id INTEGER NOT NULL, "
            "seq_region_start INTEGER NOT NULL, seq_region_end INTEGER NOT NULL, "
            "seq_region_strand INTEGER NOT NULL, repeat_start INTEGER NOT NULL, repeat_end INTEGER NOT NULL, "
            "repeat_consensus_id INTEGER NOT NULL, analysis_id INTEGER NOT NULL, score FLOAT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE attrib_type ("
            "attrib_type_id INTEGER PRIMARY KEY, code VARCHAR(20) NOT NULL UNIQUE, "
            "name VARCHAR(255) NOT NULL, description TEXT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE dna_align_feature ("
            "dna_align_feature_id INTEGER PRIMARY KEY, seq_region_id INTEGER NOT NULL, "
            "seq_region_start INTEGER NOT NULL, seq_region_end INTEGER NOT NULL, "
            "seq_region_strand INTEGER NOT NULL, "
            "hit_start INTEGER NOT NULL, hit_end INTEGER NOT NULL, hit_strand INTEGER NOT NULL, "
            "hit_name VARCHAR(40) NOT NULL, analysis_id INTEGER NOT NULL, score FLOAT, evalue FLOAT, "
            "perc_ident FLOAT, cigar_line TEXT, external_db_id INTEGER, hcoverage FLOAT, "
            "align_type VARCHAR(20))"
        )
        connection.exec_driver_sql(
            "CREATE TABLE dna_align_feature_attrib ("
            "dna_align_feature_id INTEGER NOT NULL, attrib_type_id INTEGER NOT NULL, "
            "value VARCHAR(500) NOT NULL, "
            "PRIMARY KEY (dna_align_feature_id, attrib_type_id, value))"
        )
    return Session(engine)


def _feature_values() -> dict[str, object]:
    """Return one mapped repeat feature row for persistence tests."""
    return {
        "seq_region_id": 1,
        "seq_region_start": 10,
        "seq_region_end": 20,
        "seq_region_strand": 1,
        "repeat_start": 1,
        "repeat_end": 11,
        "repeat_consensus_id": 3,
        "analysis_id": 7,
        "score": 4.5,
    }


def _alignment_values() -> dict[str, object]:
    """Return one mapped dna-align feature row for persistence tests."""
    return {
        "seq_region_id": 1,
        "seq_region_start": 10,
        "seq_region_end": 20,
        "seq_region_strand": 1,
        "hit_start": 1,
        "hit_end": 11,
        "hit_strand": 1,
        "hit_name": "RF00001",
        "analysis_id": 7,
        "score": 4.5,
        "evalue": 0.01,
        "perc_ident": None,
        "cigar_line": None,
        "external_db_id": None,
        "hcoverage": None,
        "align_type": "ensembl",
    }


def test_resolve_seq_region_id_caches_unique_names(feature_session: Session) -> None:
    """Sequence-region lookup caches each unique name and rejects invalid names.

    Args:
        feature_session: SQLite session containing feature tables.

    Raises:
        ValueError: Expected for missing or ambiguous names.
    """
    seq_region_ids: dict[str, int] = {}
    feature_session.add(SeqRegion(seq_region_id=1, name="chr1", coord_system_id=1, length=100))
    feature_session.flush()
    assert core_load._resolve_seq_region_id(feature_session, "chr1", seq_region_ids) == 1
    assert core_load._resolve_seq_region_id(feature_session, "chr1", seq_region_ids) == 1
    assert seq_region_ids == {"chr1": 1}
    with pytest.raises(ValueError, match="not found"):
        core_load._resolve_seq_region_id(feature_session, "missing", seq_region_ids)

    feature_session.add(SeqRegion(seq_region_id=2, name="chr1", coord_system_id=2, length=100))
    feature_session.flush()
    with pytest.raises(ValueError, match="ambiguous"):
        core_load._resolve_seq_region_id(feature_session, "chr1", {})


def test_get_or_create_attrib_type_reuses_and_creates(feature_session: Session) -> None:
    """Attribute types are keyed by their field code.

    Args:
        feature_session: SQLite session containing feature tables.
    """
    first_id = core_load._get_or_create_attrib_type(feature_session, "biotype")
    second_id = core_load._get_or_create_attrib_type(feature_session, "biotype")
    new_id = core_load._get_or_create_attrib_type(feature_session, "isotype")

    assert second_id == first_id
    assert new_id != first_id
    assert feature_session.query(AttribType).count() == 2


def test_get_or_create_attrib_type_rejects_long_code(feature_session: Session) -> None:
    """Core attribute type codes cannot exceed their column length.

    Args:
        feature_session: SQLite session containing feature tables.

    Raises:
        ValueError: Expected for an overlong code.
    """
    with pytest.raises(ValueError, match="longer than 20"):
        core_load._get_or_create_attrib_type(feature_session, "x" * 21)


def test_insert_repeat_feature_inserts_source_row(feature_session: Session) -> None:
    """Repeat insertion persists the source row without duplicate checking.

    Args:
        feature_session: SQLite session containing feature tables.
    """
    values = _feature_values()

    core_load._insert_repeat_feature(feature_session, values)
    assert feature_session.query(RepeatFeature).count() == 1


def test_insert_ncrna_feature_persists_attributes(
    feature_session: Session,
) -> None:
    """ncRNA insertion persists the source row and its attributes.

    Args:
        feature_session: SQLite session containing feature tables.
    """
    values = _alignment_values()
    attributes = {"biotype": "ncRNA", "ncrna_tool": "cmscan"}

    core_load._insert_ncrna_feature(feature_session, values, attributes)
    assert feature_session.query(DnaAlignFeature).count() == 1
    assert feature_session.query(DNAAlignFeatureAttrib).count() == 2


def _write_document(tmp_path: Path, document: Mapping[str, object], name: str = "input.json") -> Path:
    """Write a JSON document for loader orchestration tests.

    Args:
        tmp_path: Pytest temporary directory.
        document: JSON-compatible document to serialize.
        name: Output filename.

    Returns:
        Path to the generated document.
    """
    input_path = tmp_path / name
    input_path.write_text(json.dumps(document), encoding="utf-8")
    return input_path


def _repeat_document(feature_count: int = 2) -> dict[str, JsonValue]:
    """Return a repeat document with the requested number of features.

    Args:
        feature_count: Number of repeat features to generate.

    Returns:
        JSON-compatible repeat load document.
    """
    return {
        "analysis": _analysis_payload(),
        "source": {"source_provider": "Ensembl"},
        "repeat_features": [
            {
                "seq_region": "chr1",
                "seq_region_start": index + 1,
                "seq_region_end": index + 10,
                "seq_region_strand": "+",
                "repeat_start": 1,
                "repeat_end": 10,
            }
            for index in range(feature_count)
        ],
    }


class _RecordingSession:
    """Record loader transaction operations without connecting to MySQL."""

    def __init__(self) -> None:
        """Initialize empty transaction counters."""
        self.commit_count = 0
        self.rollback_count = 0

    def commit(self) -> None:
        """Record a simulated commit."""
        self.commit_count += 1

    def rollback(self) -> None:
        """Record a simulated rollback."""
        self.rollback_count += 1


class _RecordingDatabase:
    """Provide the session_scope interface used by DBConnectionLite."""

    def __init__(self, session: _RecordingSession) -> None:
        """Store the session returned by the context manager."""
        self.session = session

    @contextmanager
    def session_scope(self) -> Generator[_RecordingSession, None, None]:
        """Yield the recording session and roll back on errors."""
        try:
            yield self.session
        except Exception:
            self.session.rollback()
            raise


def test_iter_items_reads_array_incrementally(tmp_path: Path) -> None:
    """Array iteration yields one object at a time.

    Args:
        tmp_path: Pytest temporary directory.
    """
    input_path = _write_document(tmp_path, _repeat_document(2))
    items = core_load._iter_items(input_path, "repeat_features")

    first = next(items)
    second = next(items)

    assert first["seq_region_start"] == 1
    assert second["seq_region_start"] == 2
    with pytest.raises(StopIteration):
        next(items)


@pytest.mark.parametrize(
    ("document", "expected"),
    [
        pytest.param(_repeat_document(0), "repeat", id="repeat document"),
        pytest.param(
            {
                "analysis": _analysis_payload("cmscan_logic"),
                "ncrna_tool": "cmscan",
                "source": {"source_provider": "Ensembl"},
                "ncrna_features": [],
            },
            "ncrna",
            id="ncrna document",
        ),
    ],
)
def test_document_kind_selects_feature_array(
    tmp_path: Path, document: dict[str, object], expected: str
) -> None:
    """Document kind is selected from the top-level feature-array key.

    Args:
        tmp_path: Pytest temporary directory.
        document: JSON-compatible load document.
        expected: Expected document kind.
    """
    assert core_load._document_kind(_write_document(tmp_path, document)) == expected


def test_load_features_rejects_nonpositive_batch_size(tmp_path: Path) -> None:
    """The public loader validates its batch size independently of the CLI.

    Args:
        tmp_path: Pytest temporary directory.

    Raises:
        ValueError: Expected for a non-positive batch size.
    """
    with pytest.raises(ValueError, match="batch_size must be positive"):
        core_load.load_features(
            _write_document(tmp_path, _repeat_document(0)),
            URL.create("mysql", "user", "pw", "host", 3306, "db"),
            batch_size=0,
        )


@pytest.mark.parametrize(
    ("database_url"),
    [
        pytest.param("mysql://user:pass@localhost/db", id="mysql url"),
        pytest.param("mysql+mysqlconnector://user:pass@localhost/db", id="mysqlconnector url")
    ],
)
def test_normalise_database_url(database_url):
    """Return correct URL for mysqlconector driver.
        
    Args:
        database_url: Input database URL.
    """
    url = make_url(database_url)

    result = core_load._normalise_database_url(url)

    assert result.drivername == "mysql+mysqlconnector"


def test_lookup_consensus_id_finds_and_rejects_keys() -> None:
    """Consensus index lookups return IDs and reject unknown keys."""
    with core_load._consensus_index() as index:
        index.execute(
            "INSERT INTO consensus_key(key, consensus_id, payload) VALUES (?, ?, ?)",
            ("a" * 64, 12, "{}"),
        )
        index.commit()

        assert core_load._lookup_consensus_id(index, "a" * 64) == 12
        with pytest.raises(ValueError, match="not found"):
            core_load._lookup_consensus_id(index, "b" * 64)


def test_load_features_validates_and_commits_batches(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """The loader validates, streams, and commits each configured batch."""
    input_path = _write_document(tmp_path, _repeat_document(2))
    session = _RecordingSession()
    database = _RecordingDatabase(session)
    validation_calls: list[tuple[object, ...]] = []
    inserted: list[dict[str, object]] = []
    consensus_calls: list[dict[str, object]] = []

    def record_validation(*args: object, **kwargs: object) -> None:
        """Record schema-validator arguments."""
        validation_calls.append((*args, kwargs))

    def record_consensus(
        _session: object, consensus: dict[str, object], _key_index: sqlite3.Connection
    ) -> int:
        """Record a simulated consensus row."""
        consensus_calls.append(dict(consensus))
        return 11

    monkeypatch.setattr(core_load, "schema_validator", record_validation)
    monkeypatch.setattr(core_load, "DBConnectionLite", lambda url, reflect=False: database)
    monkeypatch.setattr(core_load, "_get_or_create_analysis", lambda session, analysis: 7)
    monkeypatch.setattr(core_load, "_get_or_create_consensus", record_consensus)
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda session, name, seq_region_ids: 1)
    monkeypatch.setattr(
        core_load,
        "_insert_repeat_feature",
        lambda session, values: inserted.append(dict(values)),
    )

    core_load.load_features(input_path, URL.create("mysql", "user", "pw", "host", 3306, "db"), batch_size=1)

    assert validation_calls == [(input_path, "load_features", {})]
    assert len(inserted) == 2
    assert len(consensus_calls) == 1
    assert session.commit_count == 2
    assert session.rollback_count == 0


def test_load_features_can_skip_validation(tmp_path: Path) -> None:
    """skip-validation bypasses schema validation while retaining feature loading."""
    input_path = _write_document(tmp_path, _repeat_document(0))
    session = _RecordingSession()
    database = _RecordingDatabase(session)
    validator = pytest.MonkeyPatch()
    validator.setattr(core_load, "schema_validator", lambda *args, **kwargs: pytest.fail("validated"))
    validator.setattr(core_load, "DBConnectionLite", lambda url, reflect=False: database)
    validator.setattr(core_load, "_get_or_create_analysis", lambda session, analysis: 7)
    validator.setattr(core_load, "_resolve_seq_region_id", lambda session, name, seq_region_ids: 1)
    try:
        core_load.load_features(
            input_path,
            URL.create("mysql", "user", "pw", "host", 3306, "db"),
            skip_validation=True,
        )
    finally:
        validator.undo()

    assert session.commit_count == 1


def test_load_features_rolls_back_active_batch_on_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Feature errors propagate and roll back the active session batch."""
    input_path = _write_document(tmp_path, _repeat_document(1))
    session = _RecordingSession()
    database = _RecordingDatabase(session)
    monkeypatch.setattr(core_load, "DBConnectionLite", lambda url, reflect=False: database)
    monkeypatch.setattr(core_load, "_get_or_create_analysis", lambda session, analysis: 7)
    monkeypatch.setattr(core_load, "_get_or_create_consensus", lambda session, consensus, key_index: 11)
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda session, name, seq_region_ids: 1)
    monkeypatch.setattr(
        core_load,
        "_insert_repeat_feature",
        lambda session, values: (_ for _ in ()).throw(ValueError("insert failed")),
    )

    with pytest.raises(ValueError, match="insert failed"):
        core_load.load_features(
            input_path,
            URL.create("mysql", "user", "pw", "host", 3306, "db"),
            skip_validation=True,
        )

    assert session.rollback_count == 1


def test_load_features_loads_ncrna_features(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """The orchestration layer dispatches ncRNA rows and commits the final batch."""
    document = {
        "analysis": _analysis_payload("ncrna_logic"),
        "ncrna_tool": "cmscan",
        "source": {"source_provider": "Ensembl"},
        "ncrna_features": [
            {
                "seq_region": "chr1",
                "seq_region_start": 1,
                "seq_region_end": 10,
                "seq_region_strand": "+",
                "biotype": "ncRNA",
                "target_name": "RF00001",
                "is_significant": True,
            }
        ],
    }
    input_path = _write_document(tmp_path, document)
    session = _RecordingSession()
    inserted: list[tuple[dict[str, object], dict[str, str]]] = []
    monkeypatch.setattr(core_load, "DBConnectionLite", lambda url, reflect=False: _RecordingDatabase(session))
    monkeypatch.setattr(core_load, "_get_or_create_analysis", lambda session, analysis: 7)
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda session, name, seq_region_ids: 1)
    monkeypatch.setattr(
        core_load,
        "_insert_ncrna_feature",
        lambda session, values, attributes: inserted.append((dict(values), dict(attributes))),
    )

    core_load.load_features(
        input_path,
        URL.create("mysql", "user", "pw", "host", 3306, "db"),
        skip_validation=True,
    )

    assert len(inserted) == 1
    assert inserted[0][0]["hit_name"] == "RF00001"
    assert inserted[0][1]["ncrna_tool"] == "cmscan"
    assert session.commit_count == 1


def test_load_features_loads_explicit_consensus_reference(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Repeat features resolve explicit consensus keys from the index."""
    key = "a" * 64
    document = _repeat_document(1)
    document["repeat_consensus"] = [_consensus_payload(key)]
    repeat_features = cast(list[dict[str, JsonValue]], document["repeat_features"])
    repeat_features[0]["repeat_consensus"] = key
    input_path = _write_document(tmp_path, document)
    session = _RecordingSession()
    consensus_calls: list[dict[str, object]] = []
    inserted: list[dict[str, object]] = []

    def record_consensus(
        _session: object, consensus: dict[str, object], key_index: sqlite3.Connection
    ) -> int:
        """Record and index a simulated consensus row."""
        consensus_calls.append(dict(consensus))
        key_index.execute(
            "INSERT INTO consensus_key(key, consensus_id, payload) VALUES (?, ?, ?)",
            (str(consensus["repeat_consensus_key"]), 11, "{}"),
        )
        key_index.commit()
        return 11

    monkeypatch.setattr(core_load, "DBConnectionLite", lambda url, reflect=False: _RecordingDatabase(session))
    monkeypatch.setattr(core_load, "_get_or_create_analysis", lambda session, analysis: 7)
    monkeypatch.setattr(core_load, "_get_or_create_consensus", record_consensus)
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda session, name, seq_region_ids: 1)
    monkeypatch.setattr(
        core_load,
        "_insert_repeat_feature",
        lambda session, values: inserted.append(dict(values)),
    )

    core_load.load_features(
        input_path,
        URL.create("mysql", "user", "pw", "host", 3306, "db"),
        skip_validation=True,
    )

    assert len(consensus_calls) == 1
    assert inserted[0]["repeat_consensus_id"] == 11


def test_main_logs_safe_url_and_forwards_options(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """The CLI hides credentials and forwards parsed options to the loader."""
    url = URL.create("mysql", "loader", "secret", "db.example", 3306, "species_core")
    args = SimpleNamespace(input_file=Path("features.json"), db_url=url, skip_validation=False, batch_size=25)
    calls: list[tuple[object, ...]] = []
    monkeypatch.setattr(core_load, "parse_args", lambda: args)
    monkeypatch.setattr(core_load, "init_logging_with_args", lambda parsed: None)
    monkeypatch.setattr(
        core_load,
        "load_features",
        lambda input_path, database_url, **kwargs: calls.append((input_path, database_url, kwargs)),
    )

    with caplog.at_level("INFO"):
        core_load.main()

    assert "secret" not in caplog.text
    assert calls == [(Path("features.json"), url, {"skip_validation": False, "batch_size": 25})]
