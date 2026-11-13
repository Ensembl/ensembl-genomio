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
"""Tests for loading GenomIO feature documents into core databases."""

# Loader behavior is intentionally tested through its private mapping helpers.
# pylint: disable=protected-access,redefined-outer-name,too-many-lines,too-many-positional-arguments

from collections.abc import Generator
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pytest
from sqlalchemy import create_engine, select
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from ensembl.core.models import (
    Analysis,
    AnalysisDescription,
    AttribType,
    DNAAlignFeatureAttrib,
    DnaAlignFeature,
    ExternalDb,
    RepeatConsensus,
    RepeatFeature,
    SeqRegion,
)
from ensembl.io.genomio.features.load import core_load

FIXTURES = Path(__file__).parent / "test_load"
CORE_URL = URL.create("mysql", "loader", "secret", "core.example", 3306, "species_core")
PROD_URL = URL.create("mysql", "reader", "secret", "prod.example", 3306, "production")


@pytest.fixture
def session() -> Generator[Session, None, None]:
    """Create the SQLite tables used by loader persistence tests.

    Yields:
        An open session backed by an in-memory SQLite database.
    """
    engine = create_engine("sqlite:///:memory:")
    with engine.begin() as connection:
        connection.exec_driver_sql(
            "CREATE TABLE analysis (analysis_id INTEGER PRIMARY KEY, created DATETIME, "
            "logic_name VARCHAR(128) NOT NULL UNIQUE, db VARCHAR(120), db_version VARCHAR(40), "
            "db_file VARCHAR(120), program VARCHAR(80), program_version VARCHAR(40), "
            "program_file VARCHAR(80), parameters TEXT, module VARCHAR(80), module_version VARCHAR(40), "
            "gff_source VARCHAR(40), gff_feature VARCHAR(40))"
        )
        connection.exec_driver_sql(
            "CREATE TABLE analysis_description (analysis_id INTEGER PRIMARY KEY, description TEXT, "
            "display_label VARCHAR(255) NOT NULL, displayable BOOLEAN DEFAULT 1, web_data TEXT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE repeat_consensus (repeat_consensus_id INTEGER PRIMARY KEY, "
            "repeat_name VARCHAR(255) NOT NULL, repeat_class VARCHAR(100) NOT NULL, "
            "repeat_type VARCHAR(40) NOT NULL, repeat_consensus TEXT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE seq_region (seq_region_id INTEGER PRIMARY KEY, name VARCHAR(255) NOT NULL, "
            "coord_system_id INTEGER NOT NULL, length INTEGER NOT NULL)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE external_db (external_db_id INTEGER PRIMARY KEY, db_name VARCHAR(100) NOT NULL, "
            "db_release VARCHAR(255), status VARCHAR(20) NOT NULL, priority INTEGER NOT NULL, "
            "db_display_name VARCHAR(255), type VARCHAR(30) NOT NULL, secondary_db_name VARCHAR(255), "
            "secondary_db_table VARCHAR(255), description TEXT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE attrib_type (attrib_type_id INTEGER PRIMARY KEY, code VARCHAR(20) NOT NULL, "
            "name VARCHAR(255) NOT NULL, description TEXT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE repeat_feature (repeat_feature_id INTEGER PRIMARY KEY, "
            "seq_region_id INTEGER NOT NULL, "
            "seq_region_start INTEGER NOT NULL, seq_region_end INTEGER NOT NULL, "
            "seq_region_strand INTEGER NOT NULL, repeat_start INTEGER NOT NULL, repeat_end INTEGER NOT NULL, "
            "repeat_consensus_id INTEGER NOT NULL, analysis_id INTEGER NOT NULL, score FLOAT)"
        )
        connection.exec_driver_sql(
            "CREATE TABLE dna_align_feature (dna_align_feature_id INTEGER PRIMARY KEY, "
            "seq_region_id INTEGER NOT NULL, seq_region_start INTEGER NOT NULL, "
            "seq_region_end INTEGER NOT NULL, "
            "seq_region_strand INTEGER NOT NULL, hit_start INTEGER NOT NULL, hit_end INTEGER NOT NULL, "
            "hit_strand INTEGER NOT NULL, hit_name VARCHAR(40) NOT NULL, analysis_id INTEGER NOT NULL, "
            "score FLOAT, evalue FLOAT, perc_ident FLOAT, cigar_line TEXT, external_db_id INTEGER, "
            "hcoverage FLOAT, align_type VARCHAR(20))"
        )
        connection.exec_driver_sql(
            "CREATE TABLE dna_align_feature_attrib (dna_align_feature_id INTEGER NOT NULL, "
            "attrib_type_id INTEGER NOT NULL, value VARCHAR(500) NOT NULL, "
            "PRIMARY KEY (dna_align_feature_id, attrib_type_id, value))"
        )
    with Session(engine) as db_session:
        yield db_session
    engine.dispose()


class RecordingSession:
    """Record commit and rollback calls made by the loader."""

    def __init__(self) -> None:
        """Initialize empty transaction counters."""
        self.commits = 0
        self.rollbacks = 0

    def commit(self) -> None:
        """Record a commit request."""
        self.commits += 1

    def rollback(self) -> None:
        """Record a rollback request."""
        self.rollbacks += 1


class RecordingDatabase:
    """Provide a ``session_scope`` context manager for loader orchestration tests."""

    def __init__(self, db_session: RecordingSession) -> None:
        """Store the session returned by the context manager.

        Args:
            db_session: Recording session to yield.
        """
        self.db_session = db_session

    @contextmanager
    def session_scope(self) -> Generator[RecordingSession, None, None]:
        """Yield the recording session and record rollback on failure.

        Yields:
            The session used by one loader invocation.
        """
        try:
            yield self.db_session
        except Exception:
            self.db_session.rollback()
            raise


@pytest.mark.parametrize(
    ("strand", "expected"),
    [("+", 1), ("-", -1), (".", 0)],
)
def test_normalise_strand(strand: str, expected: int) -> None:
    """Map GenomIO strand symbols to Ensembl integer values.

    Args:
        strand: Input strand symbol.
        expected: Expected integer strand.
    """
    assert core_load._normalise_strand(strand) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("2026-02-18T00:00:00Z", datetime(2026, 2, 18)),
        ("2026-02-18T01:00:00+01:00", datetime(2026, 2, 18)),
        ("2026-02-18T00:00:00", datetime(2026, 2, 18)),
    ],
)
def test_parse_run_date_normalises_to_utc(value: str, expected: datetime) -> None:
    """Parse ISO dates as naive UTC datetimes.

    Args:
        value: Input ISO date or datetime.
        expected: Expected naive UTC datetime.
    """
    assert core_load._parse_run_date(value) == expected


@pytest.mark.parametrize(
    ("tool", "feature", "name", "hit_start", "hit_end", "hit_strand", "cigar"),
    [
        (
            "cmscan",
            {
                "target_name": "RF00001",
                "model_start": 7,
                "model_end": 2,
                "model_strand": "-",
                "is_significant": True,
            },
            "RF00001",
            7,
            2,
            -1,
            None,
        ),
        ("cmscan", {"target_name": "RF00001", "is_significant": False}, "RF00001", 1, 21, 1, None),
        ("trnascan", {"isotype": "Phe", "anticodon": "GAA"}, "Phe", 1, 21, 1, "21M"),
        ("mirbase", {"mirna_name": "ath-miR156a"}, "ath-miR156a", 1, 21, 1, None),
        ("other", {"display_label": "label", "isotype": "X"}, "label", 1, 21, 1, None),
        ("other", {"isotype": "X"}, "X", 1, 21, 1, None),
    ],
)
def test_ncrna_row_maps_tool_specific_alignment_fields(
    monkeypatch: pytest.MonkeyPatch,
    tool: str,
    feature: dict[str, core_load.JsonValue],
    name: str,
    hit_start: int,
    hit_end: int,
    hit_strand: int,
    cigar: str | None,
) -> None:
    """Map ncRNA coordinates, names, scores, and tool-specific fields.

    Args:
        tool: Producer name.
        feature: Tool-specific input fields.
        name: Expected alignment hit name.
        hit_start: Expected hit start coordinate.
        hit_end: Expected hit end coordinate.
        hit_strand: Expected hit strand.
        cigar: Expected CIGAR string.
    """
    monkeypatch.setattr(core_load, "_get_ncrna_attributes", lambda *_: {})
    common: dict[str, core_load.JsonValue] = {
        "seq_region": "chr1",
        "seq_region_start": 10,
        "seq_region_end": 30,
        "seq_region_strand": "+",
        "score": 12,
        "evalue": 0.01,
        "biotype": "ncRNA",
    }
    common.update(feature)
    row, _ = core_load._ncrna_row(
        None,
        common,
        tool,
        seq_region_id=4,
        analysis_id=8,
        external_db_id=9,
        caches=core_load.FeatureLoadCaches(),
    )
    assert row == {
        "seq_region_id": 4,
        "seq_region_start": 10,
        "seq_region_end": 30,
        "seq_region_strand": 1,
        "hit_start": hit_start,
        "hit_end": hit_end,
        "hit_strand": hit_strand,
        "hit_name": name,
        "analysis_id": 8,
        "score": 12.0,
        "evalue": 0.01,
        "perc_ident": None,
        "cigar_line": cigar,
        "external_db_id": 9,
        "hcoverage": None,
        "align_type": "ensembl",
    }


@pytest.mark.parametrize(
    ("feature", "tool"),
    [
        ({"target_name": "x" * 41, "is_significant": True}, "cmscan"),
        ({"mirna_name": "x" * 41}, "mirbase"),
        ({"isotype": "x" * 41}, "trnascan"),
    ],
)
def test_ncrna_row_rejects_overlong_hit_names(
    monkeypatch: pytest.MonkeyPatch, feature: dict[str, core_load.JsonValue], tool: str
) -> None:
    """Reject hit names longer than the core column allows.

    Args:
        feature: Tool-specific input fields.
        tool: Producer name.

    Raises:
        ValueError: If the selected hit name exceeds 40 characters.
    """
    base: dict[str, core_load.JsonValue] = {
        "seq_region_start": 1,
        "seq_region_end": 2,
        "seq_region_strand": ".",
        "biotype": "ncRNA",
    }
    monkeypatch.setattr(core_load, "_get_ncrna_attributes", lambda *_: {})
    base.update(feature)
    with pytest.raises(ValueError, match="at most 40"):
        core_load._ncrna_row(
            None,
            base,
            tool,
            seq_region_id=1,
            analysis_id=1,
            external_db_id=None,
            caches=core_load.FeatureLoadCaches(),
        )


@pytest.mark.parametrize(
    ("tool", "feature", "codes_values"),
    [
        (
            "cmscan",
            {
                "target_accession": "RF1",
                "biotype": "ncRNA",
                "description": "desc",
                "bias": 0.5,
                "gc": 0.4,
                "truncation": "5prime",
                "is_significant": True,
            },
            [
                ("rfam_accession", "RF1"),
                ("rna_gene_biotype", "ncRNA"),
                ("description", "desc"),
                ("cmscan_bias", "0.5"),
                ("cmscan_gc", "0.4"),
                ("cmscan_truncated", "5prime"),
                ("cmscan_significant", "1"),
            ],
        ),
        (
            "cmscan",
            {"biotype": "ncRNA", "is_significant": False},
            [("rna_gene_biotype", "ncRNA"), ("cmscan_significant", "0")],
        ),
        (
            "trnascan",
            {"biotype": "tRNA", "isotype": "Phe", "anticodon": "GAA"},
            [("rna_gene_biotype", "tRNA"), ("description", "tRNA-Phe for anticodon GAA")],
        ),
        ("mirbase", {"biotype": "miRNA"}, []),
    ],
)
def test_get_ncrna_attributes_resolves_only_applicable_fields(
    session: Session,
    tool: str,
    feature: dict[str, core_load.JsonValue],
    codes_values: list[tuple[str, str]],
) -> None:
    """Resolve each applicable attribute code and serialize its value.

    Args:
        session: SQLite session containing attribute type tables.
        tool: Producer name.
        feature: Input feature fields.
        codes_values: Expected attribute code and value pairs.
    """
    for index, (code, _) in enumerate(codes_values, 1):
        session.add(AttribType(attrib_type_id=index, code=code, name=code))
    session.flush()
    attributes = core_load._get_ncrna_attributes(session, feature, tool, core_load.FeatureLoadCaches())
    assert attributes == {index: value for index, (_, value) in enumerate(codes_values, 1)}


def test_repeat_row_maps_core_columns() -> None:
    """Map repeat feature values into their core columns."""
    feature: dict[str, core_load.JsonValue] = {
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


def test_analysis_and_consensus_values_map_input_fields() -> None:
    """Map analysis and repeat consensus metadata to core column values."""
    analysis = {
        "run_date": "2026-02-18T00:00:00Z",
        "logic_name": "tool",
        "db": None,
        "db_version": "2",
        "program": "tool",
        "program_version": "1",
        "program_parameters": "-x",
    }
    assert core_load._analysis_values(analysis) == {
        "created": datetime(2026, 2, 18),
        "logic_name": "tool",
        "db": None,
        "db_version": "2",
        "program": "tool",
        "program_version": "1",
        "parameters": "-x",
    }
    assert core_load._consensus_values(
        {"repeat_name": "Alu", "repeat_class": "SINE", "repeat_type": "Alu"}
    ) == {
        "repeat_name": "Alu",
        "repeat_class": "SINE",
        "repeat_type": "Alu",
        "repeat_consensus": None,
    }


def test_get_or_create_analysis_creates_and_reuses_rows(
    session: Session, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Create analysis metadata once and reuse it on subsequent loads.

    Args:
        session: SQLite session containing analysis tables.
        monkeypatch: Fixture used to stub production metadata lookup.
    """
    monkeypatch.setattr(
        core_load,
        "get_production_analysis_description",
        lambda *_: {
            "description": "Description",
            "display_label": "Label",
            "displayable": True,
            "web_data": None,
        },
    )
    analysis = {"run_date": "2026-02-18", "logic_name": "tool", "program": "tool"}
    first = core_load._get_or_create_analysis(session, analysis, PROD_URL)
    second = core_load._get_or_create_analysis(session, analysis, PROD_URL)
    assert first == second
    assert session.scalar(select(Analysis).where(Analysis.analysis_id == first)).logic_name == "tool"
    assert (
        session.scalar(
            select(AnalysisDescription).where(AnalysisDescription.analysis_id == first)
        ).description
        == "Description"
    )


def test_get_or_create_repeat_consensus_creates_and_caches_row(session: Session) -> None:
    """Reuse matching consensus records and map each input key to its ID.

    Args:
        session: SQLite session containing repeat consensus table.
    """
    caches = core_load.FeatureLoadCaches()
    consensus = {
        "repeat_consensus_key": "key",
        "repeat_name": "Alu",
        "repeat_class": "SINE",
        "repeat_type": "Alu",
        "repeat_consensus": "ACGT",
    }
    first = core_load._get_or_create_repeat_consensus_id(session, consensus, caches)
    second = core_load._get_or_create_repeat_consensus_id(session, consensus, caches)
    assert first == second
    assert core_load._get_repeat_consensus_id("key", caches) == first
    assert session.query(RepeatConsensus).count() == 1


def test_get_or_create_repeat_consensus_reuses_matching_database_row(session: Session) -> None:
    """Reuse an existing consensus row for a new input key.

    Args:
        session: SQLite session containing repeat consensus table.
    """
    existing = RepeatConsensus(
        repeat_name="Alu", repeat_class="SINE", repeat_type="Alu", repeat_consensus="ACGT"
    )
    session.add(existing)
    session.flush()
    caches = core_load.FeatureLoadCaches()
    consensus = {
        "repeat_consensus_key": "different-key",
        "repeat_name": "Alu",
        "repeat_class": "SINE",
        "repeat_type": "Alu",
        "repeat_consensus": "ACGT",
    }
    assert (
        core_load._get_or_create_repeat_consensus_id(session, consensus, caches)
        == existing.repeat_consensus_id
    )
    assert caches.consensus_ids == {"different-key": existing.repeat_consensus_id}


def test_get_repeat_consensus_id_rejects_unknown_key() -> None:
    """Reject a feature reference whose consensus was not loaded.

    Raises:
        ValueError: If the input key is absent from the cache.
    """
    with pytest.raises(ValueError, match="key not found"):
        core_load._get_repeat_consensus_id("missing", core_load.FeatureLoadCaches())


@pytest.mark.parametrize(
    ("resolver", "model", "field", "cache_field", "name", "id_field"),
    [
        (core_load._resolve_seq_region_id, SeqRegion, "name", "seq_region_ids", "chr1", "seq_region_id"),
        (
            core_load._resolve_external_db_id,
            ExternalDb,
            "db_name",
            "external_db_ids",
            "RFAM",
            "external_db_id",
        ),
        (
            core_load._resolve_attrib_type_id,
            AttribType,
            "code",
            "attrib_type_ids",
            "biotype",
            "attrib_type_id",
        ),
    ],
)
def test_lookup_resolvers_cache_unique_database_ids(
    session: Session,
    resolver: object,
    model: type,
    field: str,
    cache_field: str,
    name: str,
    id_field: str,
) -> None:
    """Resolve and cache exact sequence, external database, and attribute names.

    Args:
        session: SQLite session containing lookup tables.
        resolver: Lookup function under test.
        model: ORM model for the lookup row.
        field: Model field matched by the name.
        cache_field: Cache dictionary field for the resolved ID.
        name: Lookup name.
        id_field: Primary key field.
    """
    extra = {"coord_system_id": 1, "length": 10} if model is SeqRegion else {}
    if model is ExternalDb:
        extra = {"status": "KNOWN", "priority": 1, "db_type": "ENSEMBL"}
    if model is AttribType:
        extra = {"name": name}
    row = model(**{id_field: 1, field: name, **extra})
    session.add(row)
    session.flush()
    caches = core_load.FeatureLoadCaches()
    assert resolver(session, name, caches) == 1
    assert resolver(session, name, caches) == 1
    assert getattr(caches, cache_field) == {name: 1}


@pytest.mark.parametrize(
    ("resolver_name", "model", "field", "name", "id_field", "extra"),
    [
        (
            "_resolve_seq_region_id",
            SeqRegion,
            "name",
            "chr1",
            "seq_region_id",
            {"coord_system_id": 1, "length": 10},
        ),
        (
            "_resolve_external_db_id",
            ExternalDb,
            "db_name",
            "RFAM",
            "external_db_id",
            {"status": "KNOWN", "priority": 1, "db_type": "ENSEMBL"},
        ),
        ("_resolve_attrib_type_id", AttribType, "code", "biotype", "attrib_type_id", {"name": "Biotype"}),
    ],
)
def test_lookup_resolvers_reject_missing_and_ambiguous_names(
    session: Session,
    resolver_name: str,
    model: type,
    field: str,
    name: str,
    id_field: str,
    extra: dict[str, object],
) -> None:
    """Reject missing names and duplicate exact matches for each lookup.

    Args:
        session: SQLite session containing lookup tables.
        resolver_name: Lookup function name.
        model: ORM row model.
        field: Name field used by lookup.
        name: Lookup key.
        id_field: Primary key field.
        extra: Other required model values.

    Raises:
        ValueError: If the lookup has zero or multiple matches.
    """
    resolver = getattr(core_load, resolver_name)
    caches = core_load.FeatureLoadCaches()
    with pytest.raises(ValueError, match="not found"):
        resolver(session, name, caches)
    session.add_all([model(**{id_field: value, field: name, **extra}) for value in (1, 2)])
    session.flush()
    with pytest.raises(ValueError, match="ambiguous"):
        resolver(session, name, caches)


def test_insert_functions_persist_rows_and_attributes(session: Session) -> None:
    """Persist repeat rows, ncRNA alignments, and their attributes.

    Args:
        session: SQLite session containing feature tables.
    """
    core_load._insert_repeat_feature(
        session,
        {
            "seq_region_id": 1,
            "seq_region_start": 1,
            "seq_region_end": 10,
            "seq_region_strand": 1,
            "repeat_start": 1,
            "repeat_end": 10,
            "repeat_consensus_id": 1,
            "analysis_id": 1,
            "score": None,
        },
    )
    core_load._insert_ncrna_feature(
        session,
        {
            "seq_region_id": 1,
            "seq_region_start": 1,
            "seq_region_end": 10,
            "seq_region_strand": 1,
            "hit_start": 1,
            "hit_end": 10,
            "hit_strand": 1,
            "hit_name": "RF00001",
            "analysis_id": 1,
            "score": None,
            "evalue": None,
            "perc_ident": None,
            "cigar_line": None,
            "external_db_id": None,
            "hcoverage": None,
            "align_type": "ensembl",
        },
        {1: "ncRNA"},
    )
    assert session.query(RepeatFeature).count() == 1
    assert session.query(DnaAlignFeature).count() == 1
    assert session.query(DNAAlignFeatureAttrib).count() == 1


@pytest.mark.parametrize(
    "filename",
    [
        "load_features_repeat.json",
        "load_features_cmscan.json",
        "load_features_trnascan.json",
        "load_features_mirbase.json",
    ],
)
def test_streaming_readers_read_existing_fixture(filename: str) -> None:
    """Read top-level feature arrays and values from schema fixture files.

    Args:
        filename: Fixture filename for one supported document kind or tool.
    """
    path = FIXTURES / filename
    kind = core_load._document_kind(path)
    assert kind == ("repeat" if "repeat" in filename else "ncrna")
    if kind == "repeat":
        assert len(list(core_load._iter_items(path, "repeat_features"))) == 1
        assert len(list(core_load._iter_items(path, "repeat_consensus"))) == 1
    else:
        assert len(list(core_load._iter_items(path, "ncrna_features"))) == 1
        assert core_load._read_value(path, "ncrna_tool") in {"cmscan", "trnascan", "mirbase"}
    assert isinstance(core_load._read_value(path, "analysis"), dict)


@pytest.mark.parametrize(
    "filename", ["load_features_cmscan.json", "load_features_trnascan.json", "load_features_mirbase.json"]
)
def test_iter_ncrna_rows_resolves_tool_database_and_maps_feature(
    session: Session,
    filename: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Yield an ncRNA row with the external database selected for its tool.

    Args:
        session: SQLite session containing sequence regions and external databases.
        filename: Existing ncRNA fixture containing the input feature.
        monkeypatch: Fixture used to record exact resolver arguments.
    """
    input_path = FIXTURES / filename
    tool = str(core_load._read_value(input_path, "ncrna_tool"))
    feature = next(core_load._iter_items(input_path, "ncrna_features"))
    db_names: list[str] = []
    expected_db = {"cmscan": "RFAM", "trnascan": "TRNASCAN_SE", "mirbase": "miRBase"}[tool]
    monkeypatch.setattr(core_load, "_iter_items", lambda *_: iter([feature]))
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda *_: 3)
    monkeypatch.setattr(
        core_load, "_resolve_external_db_id", lambda _session, name, _caches: db_names.append(name) or 5
    )
    monkeypatch.setattr(
        core_load,
        "_ncrna_row",
        lambda _session, _feature, _tool, **kwargs: ({"external_db_id": kwargs["external_db_id"]}, {}),
    )
    rows = list(core_load._iter_ncrna_rows(session, input_path, tool, 7, core_load.FeatureLoadCaches()))
    assert db_names == [expected_db]
    assert rows == [({"external_db_id": 5}, {})]


def test_parse_args_accepts_core_and_production_connection_values() -> None:
    """Parse required core and production URLs and loader options."""
    args = core_load.parse_args(
        [
            "--input_file",
            str(FIXTURES / "load_features_repeat.json"),
            "--core_host",
            "core",
            "--core_port",
            "3307",
            "--core_user",
            "loader",
            "--core_password",
            "secret",
            "--core_database",
            "core_db",
            "--prod_host",
            "prod",
            "--prod_port",
            "3308",
            "--prod_user",
            "reader",
            "--prod_password",
            "secret2",
            "--prod_database",
            "prod_db",
            "--batch-size",
            "25",
            "--skip-validation",
        ]
    )
    assert args.core_url == URL.create("mysql", "loader", "secret", "core", 3307, "core_db")
    assert args.prod_url == URL.create("mysql", "reader", "secret2", "prod", 3308, "prod_db")
    assert args.batch_size == 25
    assert args.skip_validation is True


@pytest.mark.parametrize(
    "missing",
    [
        "--core_host",
        "--core_port",
        "--core_user",
        "--core_database",
        "--prod_host",
        "--prod_port",
        "--prod_user",
        "--prod_database",
    ],
)
def test_parse_args_requires_connection_parameters(missing: str) -> None:
    """Require each core and production connection option.

    Args:
        missing: Required option removed from the argument list.

    Raises:
        SystemExit: If a required connection option is missing.
    """
    argv = ["--input_file", str(FIXTURES / "load_features_repeat.json")]
    for prefix, values in (
        ("core", ("host", "port", "user", "password", "database")),
        ("prod", ("host", "port", "user", "password", "database")),
    ):
        for value in values:
            option = f"--{prefix}_{value}"
            if option != missing:
                argv.extend((option, "3306" if value == "port" else "value"))
    with pytest.raises(SystemExit):
        core_load.parse_args(argv)


def test_load_features_rejects_nonpositive_batch_size(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reject a nonpositive batch size before validation or database setup.

    Args:
        monkeypatch: Fixture used to ensure external work is not started.

    Raises:
        ValueError: If ``batch_size`` is below one.
    """
    monkeypatch.setattr(core_load, "schema_validator", lambda *_: pytest.fail("unexpected validation"))
    with pytest.raises(ValueError, match="batch_size must be positive"):
        core_load.load_features(
            FIXTURES / "load_features_repeat.json", core_url=CORE_URL, prod_url=PROD_URL, batch_size=0
        )


@pytest.mark.parametrize("filename", ["load_features_repeat.json", "load_features_cmscan.json"])
@pytest.mark.parametrize("batch_size", [1, 2])
def test_load_features_streams_batches_and_validates(
    monkeypatch: pytest.MonkeyPatch,
    filename: str,
    batch_size: int,
) -> None:
    """Validate a fixture, load mapped rows, and commit at batch boundaries.

    Args:
        monkeypatch: Fixture used to isolate database effects.
        filename: Repeat or ncRNA fixture to load.
        batch_size: Number of features committed per batch.
    """
    input_path = FIXTURES / filename
    transaction = RecordingSession()
    inserted_repeat: list[dict[str, object]] = []
    inserted_ncrna: list[tuple[dict[str, object], dict[int, str]]] = []
    validations: list[tuple[object, ...]] = []
    monkeypatch.setattr(core_load, "schema_validator", lambda *args: validations.append(args))
    monkeypatch.setattr(core_load, "DBConnection", lambda *_args, **_kwargs: RecordingDatabase(transaction))
    monkeypatch.setattr(core_load, "_get_or_create_analysis", lambda *_: 7)
    monkeypatch.setattr(
        core_load,
        "_get_or_create_repeat_consensus_id",
        lambda _s, consensus, cache: cache.consensus_ids.update({str(consensus["repeat_consensus_key"]): 11})
        or 11,
    )
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda *_: 3)
    monkeypatch.setattr(core_load, "_resolve_external_db_id", lambda *_: 5)
    monkeypatch.setattr(core_load, "_resolve_attrib_type_id", lambda *_: 1)
    monkeypatch.setattr(
        core_load, "_insert_repeat_feature", lambda _s, row: inserted_repeat.append(dict(row))
    )
    monkeypatch.setattr(
        core_load,
        "_insert_ncrna_feature",
        lambda _s, row, attrs: inserted_ncrna.append((dict(row), dict(attrs))),
    )
    core_load.load_features(input_path, core_url=CORE_URL, prod_url=PROD_URL, batch_size=batch_size)
    assert validations == [(input_path, "load_features")]
    assert len(inserted_repeat) + len(inserted_ncrna) == 1
    assert transaction.commits == 1
    assert transaction.rollbacks == 0


def test_iter_ncrna_rows_rejects_unrecognised_tool_through_missing_database_id(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The current loader cannot construct an alignment for an unknown tool.

    Args:
        monkeypatch: Fixture used to provide a minimal input feature.

    Raises:
        UnboundLocalError: The current implementation references an unset external database ID.
    """
    monkeypatch.setattr(core_load, "_iter_items", lambda *_: iter([{"seq_region": "chr1"}]))
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda *_: 1)
    with pytest.raises(UnboundLocalError):
        list(
            core_load._iter_ncrna_rows(
                None,
                FIXTURES / "load_features_cmscan.json",
                "other",
                1,
                core_load.FeatureLoadCaches(),
            )
        )


def test_load_features_skips_validation_and_commits_empty_document(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Skip schema validation and commit when the repeat feature array is empty.

    Args:
        monkeypatch: Fixture used to isolate database effects.
    """
    transaction = RecordingSession()
    monkeypatch.setattr(core_load, "schema_validator", lambda *_: pytest.fail("unexpected validation"))
    monkeypatch.setattr(core_load, "DBConnection", lambda *_args, **_kwargs: RecordingDatabase(transaction))
    monkeypatch.setattr(core_load, "_get_or_create_analysis", lambda *_: 7)
    core_load.load_features(
        FIXTURES / "load_features_unknown_property.json",
        core_url=CORE_URL,
        prod_url=PROD_URL,
        skip_validation=True,
    )
    assert transaction.commits == 1


def test_load_features_rolls_back_on_insert_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Roll back the active transaction and propagate feature insertion errors.

    Args:
        monkeypatch: Fixture used to isolate database effects.

    Raises:
        RuntimeError: If insertion of a repeat feature fails.
    """
    transaction = RecordingSession()
    monkeypatch.setattr(core_load, "schema_validator", lambda *_: None)
    monkeypatch.setattr(core_load, "DBConnection", lambda *_args, **_kwargs: RecordingDatabase(transaction))
    monkeypatch.setattr(core_load, "_get_or_create_analysis", lambda *_: 7)
    monkeypatch.setattr(
        core_load,
        "_get_or_create_repeat_consensus_id",
        lambda _s, consensus, cache: cache.consensus_ids.update({str(consensus["repeat_consensus_key"]): 11})
        or 11,
    )
    monkeypatch.setattr(core_load, "_resolve_seq_region_id", lambda *_: 3)
    monkeypatch.setattr(
        core_load, "_insert_repeat_feature", lambda *_: (_ for _ in ()).throw(RuntimeError("insert failed"))
    )
    with pytest.raises(RuntimeError, match="insert failed"):
        core_load.load_features(FIXTURES / "load_features_repeat.json", core_url=CORE_URL, prod_url=PROD_URL)
    assert transaction.rollbacks == 1


def test_main_uses_driver_urls_and_logs_without_password(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Forward CLI options with MySQL connector URLs and hide credentials in logs.

    Args:
        monkeypatch: Fixture used to stub CLI and database operations.
        caplog: Fixture capturing log output.
    """
    args = SimpleNamespace(
        input_file=FIXTURES / "load_features_repeat.json",
        core_url=CORE_URL,
        prod_url=PROD_URL,
        skip_validation=True,
        batch_size=25,
    )
    received: list[tuple[object, ...]] = []
    monkeypatch.setattr(core_load, "parse_args", lambda _argv=None: args)
    monkeypatch.setattr(core_load, "load_features", lambda *a, **kw: received.append((a, kw)))
    with caplog.at_level("INFO"):
        core_load.main([])
    assert "secret" not in caplog.text
    positional, keyword = received[0]
    assert positional == (args.input_file,)
    assert keyword == {
        "core_url": CORE_URL.set(drivername="mysql+mysqlconnector"),
        "prod_url": PROD_URL.set(drivername="mysql+mysqlconnector"),
        "skip_validation": True,
        "batch_size": 25,
    }


def test_main_logs_and_reraises_loader_error(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Log loader failures and re-raise the original exception.

    Args:
        monkeypatch: Fixture used to stub CLI parsing and loader behavior.
        caplog: Fixture capturing error log output.

    Raises:
        RuntimeError: Propagated from the loader.
    """
    monkeypatch.setattr(
        core_load,
        "parse_args",
        lambda _argv=None: SimpleNamespace(
            input_file=Path("features.json"),
            core_url=CORE_URL,
            prod_url=PROD_URL,
            skip_validation=False,
            batch_size=1000,
        ),
    )
    monkeypatch.setattr(
        core_load, "load_features", lambda *_a, **_kw: (_ for _ in ()).throw(RuntimeError("failed"))
    )
    with caplog.at_level("ERROR"), pytest.raises(RuntimeError, match="failed"):
        core_load.main([])
    assert "Error occurred while loading features" in caplog.text
