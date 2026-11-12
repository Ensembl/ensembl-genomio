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
"""Load GenomIO repeat and ncRNA feature JSON into an Ensembl core database."""

import argparse
from collections.abc import Generator, Iterator, Mapping
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import logging
from os import name
from pathlib import Path
import sqlite3
import tempfile
from typing import TypedDict, cast


import ijson
from python.ensembl.io.genomio.database.analysis import get_production_analysis_description
from sqlalchemy import select
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

import ensembl.io.genomio
from ensembl.io.genomio.database import get_production_analysis_description
from ensembl.io.genomio.schemas.json.validate import schema_validator

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.database import DBConnection
from ensembl.utils.logging import init_logging_with_args

__all__ = ["load_features", "main", "parse_args"]

JsonValue = bool | int | float | str | list["JsonValue"] | dict[str, "JsonValue"] | None

MAX_ATTRIB_TYPE_CODE_LENGTH = 20
MAX_NCRNA_HIT_NAME_LENGTH = 40


class DnaAlignFeatureRow(TypedDict):
    """Values stored in a dna_align_feature row."""

    seq_region_id: int
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: int
    hit_start: int
    hit_end: int
    hit_strand: int
    hit_name: str
    analysis_id: int
    score: float | None
    evalue: float | None
    perc_ident: float | None
    cigar_line: str | None
    external_db_id: int | None
    hcoverage: float | None
    align_type: str


class RepeatFeatureRow(TypedDict):
    """Values stored in a repeat_feature row."""

    seq_region_id: int
    seq_region_start: int
    seq_region_end: int
    seq_region_strand: int
    repeat_start: int
    repeat_end: int
    repeat_consensus_id: int
    analysis_id: int
    score: float | None


def _normalise_strand(strand: str) -> int:
    """Convert a GenomIO strand into the Ensembl integer representation."""
    return {"+": 1, "-": -1, ".": 0}[strand]


def _attribute_scalar_to_text(value: object) -> str:
    """Convert a schema scalar to the text required by attribute values."""
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def _ncrna_hit_name(feature: Mapping[str, JsonValue], tool: str) -> str:
    """Select and validate the Ensembl hit name for an ncRNA feature.

    Args:
        feature: One validated ncRNA feature object.
        tool: Producing tool that determines the source field.

    Returns:
        A hit name no longer than the core column limit.

    Raises:
        ValueError: If the selected name exceeds the core column length.

    """
    if tool == "cmscan":
        name = feature["target_name"]
    elif tool == "mirbase":
        name = feature["mirna_name"]
    else:
        name = feature.get("display_label") or feature["isotype"]
    name = cast("str", name)
    if len(name) > MAX_NCRNA_HIT_NAME_LENGTH:
        raise ValueError(f"ncRNA hit_name must be at most {MAX_NCRNA_HIT_NAME_LENGTH} characters: {name!r}")
    return name


def _ncrna_row(
    feature: Mapping[str, JsonValue],
    tool: str,
    seq_region_id: int,
    analysis_id: int,
    external_db_id: int | None,
) -> tuple[DnaAlignFeatureRow, dict[str, str]]:
    """Map one schema ncRNA object to a row and its feature attributes.

    Args:
        feature: One validated ncRNA feature object.
        tool: Producing tool, used to select hit-name and coordinate rules.
        seq_region_id: Database identifier resolved from ``feature["seq_region"]``.
        analysis_id: Database identifier for the owning analysis.
        external_db_id: Database identifier for the external database, or None if not applicable.

    Returns:
        The ``dna_align_feature`` values and serialized attribute values.

    Raises:
        ValueError: If the selected hit name exceeds the core column length.

    """
    genomic_length = cast("int", feature["seq_region_end"]) - cast("int", feature["seq_region_start"]) + 1
    if tool == "cmscan":
        hit_start = cast("int", feature.get("model_start") or 1)
        hit_end = cast("int", feature.get("model_end") or genomic_length)
    else:
        hit_start = 1
        hit_end = genomic_length

    if tool == "trnascan":
        cigar_line = f"{genomic_length}M"
    else:
        cigar_line = None

    direct_fields = {
        "seq_region",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "score",
        "evalue",
        "target_name",
        "mirna_name",
        "model_start",
        "model_end",
    }
    attributes = {
        key: _attribute_scalar_to_text(value) for key, value in feature.items() if key not in direct_fields
    }
    attributes["ncrna_tool"] = tool
    row: DnaAlignFeatureRow = {
        "seq_region_id": seq_region_id,
        "seq_region_start": cast("int", feature["seq_region_start"]),
        "seq_region_end": cast("int", feature["seq_region_end"]),
        "seq_region_strand": _normalise_strand(str(feature["seq_region_strand"])),
        "hit_start": hit_start,
        "hit_end": hit_end,
        "hit_strand": 1,
        "hit_name": _ncrna_hit_name(feature, tool),
        "analysis_id": analysis_id,
        "score": float(cast("int | float", feature["score"])) if "score" in feature else None,
        "evalue": float(cast("int | float", feature["evalue"])) if "evalue" in feature else None,
        "perc_ident": None,
        "cigar_line": cigar_line,
        "external_db_id": external_db_id,
        "hcoverage": None,
        "align_type": "ensembl",
    }
    return row, attributes


def _parse_run_date(value: str) -> datetime:
    """Parse an ISO run date and normalize it to UTC-naive form.

    Args:
        value: ISO-8601 date or datetime string.

    Returns:
        A naive UTC ``datetime`` suitable for the core ``analysis.created`` column.

    Raises:
        ValueError: If ``value`` is not a valid ISO-8601 datetime.

    """
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        return parsed
    return parsed.astimezone(timezone.utc).replace(tzinfo=None)


def _repeat_row(
    feature: Mapping[str, JsonValue],
    analysis_id: int,
    seq_region_id: int,
    consensus_id: int,
) -> RepeatFeatureRow:
    """Map one schema repeat object to a ``repeat_feature`` row.

    Args:
        feature: One validated repeat feature object.
        analysis_id: Database identifier for the owning analysis.
        seq_region_id: Database identifier for the feature's sequence region.
        consensus_id: Database identifier for the repeat consensus.

    Returns:
        Values for one ``repeat_feature`` insert.

    """
    row: RepeatFeatureRow = {
        "seq_region_id": seq_region_id,
        "seq_region_start": cast("int", feature["seq_region_start"]),
        "seq_region_end": cast("int", feature["seq_region_end"]),
        "seq_region_strand": _normalise_strand(str(feature["seq_region_strand"])),
        "repeat_start": cast("int", feature["repeat_start"]),
        "repeat_end": cast("int", feature["repeat_end"]),
        "repeat_consensus_id": consensus_id,
        "analysis_id": analysis_id,
        "score": float(cast("int | float", feature["score"])) if "score" in feature else None,
    }

    return row


def _analysis_values(analysis: Mapping[str, JsonValue]) -> dict[str, object]:
    """Return the values persisted in the ``analysis`` table.

    Args:
        analysis: Validated analysis metadata from the input document.

    Returns:
        Core ``analysis`` column values.

    """
    return {
        "created": _parse_run_date(str(analysis["run_date"])),
        "logic_name": str(analysis["logic_name"]),
        "db": analysis.get("db"),
        "db_version": analysis.get("db_version"),
        "program": analysis.get("program"),
        "program_version": analysis.get("program_version"),
        "parameters": analysis.get("program_parameters"),
    }


def _get_or_create_analysis(session: Session, analysis: Mapping[str, JsonValue], prod_url: URL) -> int:
    """Create or verify the analysis and its description.

    Args:
        session: Active core database session.
        analysis: Validated analysis metadata from the input document.
        prod_url: URL for the production database.

    Returns:
        The existing or newly created analysis ID.

    """
    values = _analysis_values(analysis)
    existing = session.scalar(select(Analysis).where(Analysis.logic_name == values["logic_name"]))
    if existing is not None:
        return int(existing.analysis_id)

    prod_analysis_description_values = get_production_analysis_description(
        str(values["logic_name"]), prod_url
    )
    created = Analysis(**values)
    session.add(created)
    session.flush()
    session.add(
        AnalysisDescription(
            analysis_id=created.analysis_id,
            description=prod_analysis_description_values["description"],
            display_label=prod_analysis_description_values["display_label"],
            displayable=prod_analysis_description_values["displayable"],
            web_data=prod_analysis_description_values["web_data"],
        )
    )
    session.flush()
    return int(created.analysis_id)


@contextmanager
def _consensus_index() -> Generator[sqlite3.Connection, None, None]:
    """Create and remove a temporary disk-backed consensus key index.

    Yields:
        A SQLite connection containing the input-key to core-ID index.

    """
    with tempfile.TemporaryDirectory(prefix="genomio-consensus-") as directory:
        connection = sqlite3.connect(Path(directory) / "index.sqlite3")
        connection.execute(
            "CREATE TABLE consensus_key "
            "(key TEXT PRIMARY KEY, consensus_id INTEGER NOT NULL, payload TEXT NOT NULL)"
        )
        try:
            yield connection
        finally:
            connection.close()


def _consensus_values(consensus: Mapping[str, JsonValue]) -> dict[str, object]:
    """Return values persisted in ``repeat_consensus``.

    Args:
        consensus: Validated repeat consensus object.

    Returns:
        Core ``repeat_consensus`` column values.

    """
    return {
        "repeat_name": consensus["repeat_name"],
        "repeat_class": consensus["repeat_class"],
        "repeat_type": consensus["repeat_type"],
        "repeat_consensus": consensus.get("repeat_consensus"),
    }


def _get_or_create_consensus(
    session: Session,
    consensus: Mapping[str, JsonValue],
    key_index: sqlite3.Connection,
) -> int:
    """Create or reuse a repeat consensus and index its JSON key.

    Args:
        session: Active core database session.
        consensus: Repeat consensus object containing its input key.
        key_index: Temporary disk-backed input-key index.

    Returns:
        The existing or newly created repeat consensus ID.

    Raises:
        ValueError: If one input key is associated with conflicting data.

    """
    key = str(consensus["repeat_consensus_key"])
    values = _consensus_values(consensus)
    payload = json.dumps(values, sort_keys=True, separators=(",", ":"))
    indexed = key_index.execute(
        "SELECT consensus_id, payload FROM consensus_key WHERE key = ?", (key,)
    ).fetchone()
    if indexed is not None:
        if indexed[1] != payload:
            raise ValueError(f"Repeat consensus key {key} has conflicting consensus data")
        return int(indexed[0])

    filters = [
        (
            getattr(RepeatConsensus, field).is_(value)
            if value is None
            else getattr(RepeatConsensus, field) == value
        )
        for field, value in values.items()
    ]
    existing = session.scalar(select(RepeatConsensus).where(*filters))
    if existing is None:
        existing = RepeatConsensus(**values)
        session.add(existing)
        session.flush()
    key_index.execute(
        "INSERT INTO consensus_key(key, consensus_id, payload) VALUES (?, ?, ?)",
        (key, int(existing.repeat_consensus_id), payload),
    )
    key_index.commit()
    return int(existing.repeat_consensus_id)


def _resolve_seq_region_id(session: Session, name: str, seq_region_ids: dict[str, int]) -> int:
    """Resolve and cache one exact sequence-region name to its database ID.

    Args:
        session: Active core database session.
        name: Sequence region name.
        seq_region_ids: Per-load sequence-region name to ID cache.

    Returns:
        The unique matching sequence-region ID.

    Raises:
        ValueError: If no unique sequence region matches ``name``.

    """
    if name in seq_region_ids:
        return seq_region_ids[name]
    rows = session.scalars(select(SeqRegion).where(SeqRegion.name == name)).all()
    if not rows:
        raise ValueError(f"Sequence region not found: {name}")
    if len(rows) > 1:
        raise ValueError(f"Sequence region name is ambiguous: {name}")
    seq_region_id = int(rows[0].seq_region_id)
    seq_region_ids[name] = seq_region_id
    return seq_region_id


def _resolve_external_db_id(session: Session, db_name: str, external_db_ids: dict[str, int]) -> int:
    """Resolve and cache one exact external database name to its database ID.

    Args:
        session: Active core database session.
        name: External database name.
        external_db_ids: Per-load external database name to ID cache.

    Returns:
        The unique matching external database ID.

    Raises:
        ValueError: If no unique external database matches ``name``.

    """
    if db_name in external_db_ids:
        return external_db_ids[db_name]
    rows = session.scalars(select(ExternalDb).where(ExternalDb.db_name == db_name)).all()
    if not rows:
        raise ValueError(f"External database not found: {db_name}")
    if len(rows) > 1:
        raise ValueError(f"External database name is ambiguous: {db_name}")
    external_db_id = int(rows[0].external_db_id)
    external_db_ids[db_name] = external_db_id
    return external_db_id


def _get_or_create_attrib_type(session: Session, code: str) -> int:
    """Return the core attribute type ID for one ncRNA field name.

    Args:
        session: Active core database session.
        code: Attribute code, limited by the core column to {MAX_ATTRIB_TYPE_CODE_LENGTH} characters.

    Returns:
        The existing or newly created attribute type ID.

    Raises:
        ValueError: If ``code`` exceeds the core column length.

    """
    if len(code) > MAX_ATTRIB_TYPE_CODE_LENGTH:
        raise ValueError(
            f"Attribute type code is longer than {MAX_ATTRIB_TYPE_CODE_LENGTH} characters: {code}"
        )
    existing = session.scalar(select(AttribType).where(AttribType.code == code))
    if existing is None:
        existing = AttribType(code=code, name=code, description="GenomIO ncRNA feature attribute")
        session.add(existing)
        session.flush()
    return int(existing.attrib_type_id)


def _insert_repeat_feature(session: Session, values: Mapping[str, object]) -> None:
    """Insert one repeat feature from the source document.

    Args:
        session: Active core database session.
        values: ``repeat_feature`` column values.

    """
    session.add(RepeatFeature(**values))
    session.flush()


def _insert_ncrna_feature(
    session: Session,
    values: Mapping[str, object],
    attributes: Mapping[str, str],
) -> None:
    """Insert one ncRNA alignment and its source attributes.

    Args:
        session: Active core database session.
        values: ``dna_align_feature`` column values.
        attributes: Attribute code/value pairs to store with the feature.

    """
    feature = DnaAlignFeature(**values)
    session.add(feature)
    session.flush()
    for code, value in attributes.items():
        session.add(
            DNAAlignFeatureAttrib(
                dna_align_feature_id=feature.dna_align_feature_id,
                attrib_type_id=_get_or_create_attrib_type(session, code),
                value=value,
            )
        )
    session.flush()


def _iter_items(input_path: Path, prefix: str) -> Iterator[dict[str, JsonValue]]:
    """Yield objects from one top-level JSON array without retaining the array.

    Args:
        input_path: Input JSON path.
        prefix: Top-level array name to stream.

    Yields:
        One decoded array item at a time.

    """
    with Path(input_path).open("rb") as input_handle:
        yield from ijson.items(input_handle, f"{prefix}.item", use_float=True)


def _read_value(input_path: Path, prefix: str) -> object:
    """Read one small top-level JSON value from a streaming parser.

    Args:
        input_path: Input JSON path.
        prefix: Top-level value name to read.

    Returns:
        The decoded value.

    """
    with Path(input_path).open("rb") as input_handle:
        return next(ijson.items(input_handle, prefix, use_float=True))


def _document_kind(input_path: Path) -> str:
    """Return whether a document contains repeat or ncRNA features.

    Args:
        input_path: Input JSON path.

    Returns:
        ``"repeat"`` or ``"ncrna"``.

    """
    feature_keys: set[str] = set()
    with Path(input_path).open("rb") as input_handle:
        for prefix, event, value in ijson.parse(input_handle, use_float=True):
            if prefix == "" and event == "map_key" and value in {"repeat_features", "ncrna_features"}:
                feature_keys.add(str(value))
    return "repeat" if "repeat_features" in feature_keys else "ncrna"


def _lookup_consensus_id(key_index: sqlite3.Connection, key: str) -> int:
    """Resolve a consensus key from the disk-backed index.

    Args:
        key_index: Temporary consensus-key SQLite connection.
        key: Input repeat consensus key.

    Returns:
        The corresponding core repeat consensus ID.

    Raises:
        ValueError: If ``key`` has not been loaded.

    """
    row = key_index.execute("SELECT consensus_id FROM consensus_key WHERE key = ?", (key,)).fetchone()
    if row is None:
        raise ValueError(f"Repeat consensus key not found: {key}")
    return int(row[0])


def _iter_repeat_rows(
    session: Session,
    input_path: Path,
    *,
    analysis_id: int,
    key_index: sqlite3.Connection,
    seq_region_ids: dict[str, int],
) -> Iterator[tuple[Mapping[str, object], None]]:
    """Yield mapped repeat rows while resolving their consensus IDs.

    Args:
        session: Active core database session.
        input_path: Input JSON path.
        analysis_id: Database identifier for the owning analysis.
        logic_name: Analysis logic name used for synthetic consensus rows.
        key_index: Temporary consensus-key SQLite connection.
        seq_region_ids: Per-load sequence-region name to ID cache.

    Yields:
        Repeat feature values paired with ``None`` attribute values.

    """
    for feature in _iter_items(input_path, "repeat_features"):
        seq_region_id = _resolve_seq_region_id(session, str(feature["seq_region"]), seq_region_ids)
        consensus_key = feature.get("repeat_consensus")
        consensus_id = _lookup_consensus_id(key_index, str(consensus_key))
        yield _repeat_row(feature, analysis_id, seq_region_id, consensus_id), None


def _iter_ncrna_rows(
    session: Session,
    input_path: Path,
    tool: str,
    analysis_id: int,
    seq_region_ids: dict[str, int],
    external_db_ids: dict[str, int],
) -> Iterator[tuple[Mapping[str, object], dict[str, str]]]:
    """Yield mapped ncRNA rows and their attribute values.

    Args:
        session: Active core database session.
        input_path: Input JSON path.
        tool: Producing ncRNA tool.
        analysis_id: Database identifier for the owning analysis.
        seq_region_ids: Per-load sequence-region name to ID cache.

    Yields:
        Alignment values paired with their attribute code/value pairs.

    """
    for feature in _iter_items(input_path, "ncrna_features"):
        seq_region_id = _resolve_seq_region_id(session, str(feature["seq_region"]), seq_region_ids)
        external_db_id = _resolve_external_db_id(session, str(feature[""]), external_db_ids)
        values, attributes = _ncrna_row(feature, tool, seq_region_id, analysis_id, external_db_id)
        yield values, attributes


def load_features(
    input_path: Path,
    *,
    core_url: URL,
    prod_url: URL,
    skip_validation: bool = False,
    batch_size: int = 1000,
) -> None:
    """Stream source-of-truth feature rows into an Ensembl core database.

    Args:
        input_path: JSON document conforming to ``load_features.json``.
        database_url: SQLAlchemy URL for the target MySQL core database.
        skip_validation: Whether to skip validation of a document already validated by the caller.
        batch_size: Number of input features between explicit commits.

    Raises:
        ValueError: If the batch size is not positive.

    """
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if not skip_validation:
        schema_validator(input_path, "load_features")

    kind = _document_kind(input_path)
    analysis = cast("Mapping[str, JsonValue]", _read_value(input_path, "analysis"))
    tool = ""
    if kind == "ncrna":
        tool = cast("str", _read_value(input_path, "ncrna_tool"))

    database = DBConnection(core_url, reflect=False)
    with _consensus_index() as key_index, database.session_scope() as session:
        seq_region_ids: dict[str, int] = {}
        analysis_id = _get_or_create_analysis(session, analysis, prod_url)
        rows: Iterator[tuple[Mapping[str, object], None | dict[str, str]]]
        if kind == "repeat":
            for consensus in _iter_items(input_path, "repeat_consensus"):
                _get_or_create_consensus(session, consensus, key_index)

        if kind == "repeat":
            rows = _iter_repeat_rows(
                session=session,
                input_path=input_path,
                analysis_id=analysis_id,
                key_index=key_index,
                seq_region_ids=seq_region_ids,
            )
        else:
            external_db_ids: dict[str, int] = {}
            rows = _iter_ncrna_rows(
                session=session,
                input_path=input_path,
                tool=tool,
                analysis_id=analysis_id,
                seq_region_ids=seq_region_ids,
                external_db_ids=external_db_ids,
            )

        feature_count = 0
        for values, attributes in rows:
            if attributes is None:
                _insert_repeat_feature(session, values)
            else:
                _insert_ncrna_feature(session, values, attributes)
            feature_count += 1
            if feature_count % batch_size == 0:
                session.commit()

        if feature_count % batch_size or feature_count == 0:
            session.commit()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments using required command-line database values.

    Args:
        arg_list: Optional argument list; defaults to the process command line.

    Returns:
        Parsed command-line options, including the constructed database URL.

    """
    parser = ArgumentParser(description="Load GenomIO features into an Ensembl core database.")
    parser.add_argument_src_path("--input_file", required=True, help="load_features JSON document.")
    parser.add_argument("--skip-validation", action="store_true", help="skip JSON schema validation.")
    parser.add_numeric_argument(
        "--batch-size",
        type=int,
        min_value=1,
        default=1000,
        help="number of features to commit per batch",
    )
    parser.add_server_arguments(
        prefix="core_",
        include_database=True,
        help="Core database connection for loading features.",
    )
    parser.add_server_arguments(
        prefix="prod_",
        include_database=True,
        help="Production database connection for retrieving analysis metadata.",
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()

    args = parser.parse_args(argv)
    init_logging_with_args(args)
    return args


def main(argv: list[str] | None = None) -> None:
    """Run the features loader CLI.

    Args:
        argv: Optional list of command-line arguments; defaults to the process command line.

    """
    args = parse_args(argv)
    logging.info("Loading features into %s", args.db_url.render_as_string(hide_password=True))

    try:
        load_features(
            args.input_file,
            core_url=args.core_url.set(drivername="mysql+mysqlconnector"),
            prod_url=args.prod_url.set(drivername="mysql+mysqlconnector"),
            skip_validation=args.skip_validation,
            batch_size=args.batch_size,
        )
    except Exception:
        logging.exception(f"Error occurred while loading features from {args.input_file}")
        raise
