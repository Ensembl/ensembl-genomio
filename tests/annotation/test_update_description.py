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
"""Unit testing of `ensembl.io.genomio.annotation.update_description` module."""
# pylint: disable=too-many-positional-arguments

from contextlib import nullcontext as no_raise
from pathlib import Path
from typing import Callable, ContextManager

import pytest
from pytest import CaptureFixture, param, raises
import sqlalchemy
from sqlalchemy import text

from ensembl.io.genomio.annotation import get_core_data, load_descriptions
from ensembl.core.models import Gene, Transcript, ObjectXref, Xref, metadata
from ensembl.utils.database import UnitTestDB


def add_gene(dialect: str, session: sqlalchemy.orm.Session, gene_data: dict[str, str]) -> None:
    """Add a gene and a transcript from a gene_data dict.

    Args:
    dialect: if "mysql", allow to skip foreign key check.
    session: SQLalchemy session.
    gene_data:
        "gene_name" -> gene.stable_id
        "gene_desc" -> gene.description
        "tr_name" -> transcript.stable_id
        "tr_desc" -> transcript.description
        "gene_xref" -> xref display_name attached to the gene
    """
    gene_name = gene_data.get("gene_name", "gene1")
    gene_description = gene_data.get("gene_desc", "")
    gene_xref = gene_data.get("gene_xref", "")
    tr_name = gene_data.get("tr_name", "tr1")
    tr_description = gene_data.get("tr_desc", "")
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=0"))
    new_transcript = Transcript(
        gene_id=1,
        transcript_id=1,
        stable_id=tr_name,
        biotype="protein_coding",
        description=tr_description,
        analysis_id=1,
        seq_region_id=1,
        seq_region_start=1,
        seq_region_end=100,
        seq_region_strand=1,
        source="test_source",
    )

    new_gene = Gene(
        gene_id=1,
        stable_id=gene_name,
        biotype="protein_coding",
        description=gene_description,
        analysis_id=1,
        seq_region_id=1,
        seq_region_start=1,
        seq_region_end=100,
        seq_region_strand=1,
        source="test_source",
        canonical_transcript_id=1,
    )
    session.add(new_gene)
    session.add(new_transcript)

    if gene_xref:
        new_xref = Xref(xref_id=1, dbprimary_acc=gene_xref, display_label=gene_xref, external_db_id=1)
        new_object_xref = ObjectXref(xref_id=1, ensembl_object_type="gene", ensembl_id=1)
        session.add(new_xref)
        session.add(new_object_xref)

    session.commit()
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=1"))


@pytest.fixture(name="annot_test_db", scope="module")
def fixture_annotation_test_db(db_factory: Callable) -> UnitTestDB:
    """Returns a test database with all core tables and basic data."""
    test_db: UnitTestDB = db_factory(src="no_src", name="load_annotation")
    test_db.dbc.create_all_tables(metadata)
    return test_db


@pytest.mark.parametrize(
    "gene_data, table, expected_ids, expectation",
    [
        param({"name": "gene1", "tr_name": "tr1"}, "gene", ["gene1"], no_raise()),
        param({"name": "gene1", "tr_name": "tr1"}, "transcript", ["tr1"], no_raise()),
        param({}, "foo", [], raises(ValueError, match="Table foo is not supported")),
    ],
)
def test_get_core_data(
    annot_test_db: UnitTestDB,
    gene_data: dict[str, str],
    table: str,
    expected_ids: list[str],
    expectation: ContextManager,
) -> None:
    """Tests the method `get_core_data()`"""
    with annot_test_db.dbc.test_session_scope() as session:
        add_gene(annot_test_db.dbc.dialect, session, gene_data)
        with expectation:
            feats = get_core_data(session, table)
            assert list(feats.keys()) == expected_ids


@pytest.mark.parametrize(
    "input_file, gene_data, do_update, expected_description",
    [
        param("gene1_desc.json", {"gene_name": "gene1"}, False, "", id="No update"),
        param("gene1_desc.json", {"gene_name": "gene1"}, True, "new_desc", id="Do update"),
    ],
)
def test_load_description_do_update(
    annot_test_db: UnitTestDB,
    data_dir: Path,
    input_file: str,
    gene_data: dict[str, str],
    do_update: bool,
    expected_description: str,
) -> None:
    """Tests the method `load_description()`"""
    with annot_test_db.dbc.test_session_scope() as session:
        add_gene(annot_test_db.dbc.dialect, session, gene_data)
        load_descriptions(session, data_dir / input_file, do_update=do_update)
        feats = get_core_data(session, "gene")
        assert len(feats) == 1
        name = gene_data["gene_name"]
        assert feats[name]
        assert feats[name][2] == expected_description


@pytest.mark.parametrize(
    "input_file, gene_data, table, expected_description",
    [
        param("gene1_nodesc.json", {"gene_name": "gene1"}, "gene", "", id="Gene: no desc -> no desc"),
        param(
            "gene1_nodesc.json",
            {"gene_name": "gene1", "gene_desc": "old_desc"},
            "gene",
            "",
            id="Gene: old desc -> deleted desc",
        ),
        param("gene1_desc.json", {"gene_name": "gene1"}, "gene", "new_desc", id="Gene: no desc -> new desc"),
        param(
            "gene1_desc.json",
            {"gene_name": "gene1", "gene_desc": "new_desc"},
            "gene",
            "new_desc",
            id="Gene: old desc -> same desc",
        ),
        param(
            "gene1_desc.json",
            {"gene_name": "gene1", "gene_desc": "old_desc"},
            "gene",
            "new_desc",
            id="Gene: old desc -> new desc",
        ),
        param(
            "gene1_nodesc.json",
            {"gene_name": "gene1", "gene_desc": "old_desc [Source: ext]"},
            "gene",
            "old_desc [Source: ext]",
            id="Gene: source desc -> no change",
        ),
        param(
            "gene1_desc.json",
            {"gene_name": "gene1", "gene_desc": "old_desc [Source: ext]"},
            "gene",
            "new_desc",
            id="Gene: source desc -> new desc",
        ),
        param(
            "gene1_syn.json",
            {"gene_name": "gene1"},
            "gene",
            "new_desc",
            id="Gene: no desc -> new_desc from syn match",
        ),
        param(
            "gene1_syn_nouse.json",
            {"gene_name": "gene1"},
            "gene",
            "",
            id="Gene: no desc -> no syn match",
        ),
        param(
            "gene1_xref.json",
            {"gene_name": "gene1"},
            "gene",
            "new_desc",
            id="Gene: no desc -> new_desc from xref match",
        ),
        param(
            "gene1_xref_nouse.json",
            {"gene_name": "gene1"},
            "gene",
            "",
            id="Gene: no desc -> no xref match",
        ),
        param("tr1_nodesc.json", {"tr_name": "tr1"}, "transcript", "", id="Tr: no desc -> no desc"),
        param("tr1_desc.json", {"tr_name": "tr1"}, "transcript", "new_desc", id="Tr: no desc -> new desc"),
        param(
            "tr1_desc.json",
            {"tr_name": "tr1", "tr_desc": "old_desc"},
            "transcript",
            "new_desc",
            id="Tr: old desc -> new desc",
        ),
    ],
)
def test_load_description(
    annot_test_db: UnitTestDB,
    data_dir: Path,
    input_file: str,
    gene_data: dict[str, str],
    table: str,
    expected_description: str,
) -> None:
    """Tests the method `load_description()`"""
    with annot_test_db.dbc.test_session_scope() as session:
        add_gene(annot_test_db.dbc.dialect, session, gene_data)
        load_descriptions(session, data_dir / input_file, do_update=True)
        feats = get_core_data(session, table)
        assert len(feats) == 1
        name = ""
        if table == "gene":
            name = gene_data["gene_name"]
        elif table == "transcript":
            name = gene_data["tr_name"]
        assert feats[name]
        assert feats[name][2] == expected_description


@pytest.mark.parametrize(
    "input_file, gene_data, expected_description, match_xrefs",
    [
        param(
            "gene1_desc.json",
            {"gene_name": "foobar", "gene_xref": "gene1"},
            "",
            False,
            id="Gene: no desc -> no changed no xref match",
        ),
        param(
            "gene1_desc.json",
            {"gene_name": "foobar", "gene_xref": "gene1"},
            "new_desc",
            True,
            id="Gene: no desc -> new_desc from xref match from json",
        ),
    ],
)
def test_load_description_match_xrefs(
    annot_test_db: UnitTestDB,
    data_dir: Path,
    input_file: str,
    gene_data: dict[str, str],
    match_xrefs: bool,
    expected_description: str,
) -> None:
    """Tests the method `load_description()` with `match_xrefs`"""
    with annot_test_db.dbc.test_session_scope() as session:
        add_gene(annot_test_db.dbc.dialect, session, gene_data)
        load_descriptions(session, data_dir / input_file, do_update=True, match_xrefs=match_xrefs)
        feats = get_core_data(session, "gene")
        name = gene_data["gene_name"]
        assert feats[name]
        assert feats[name][2] == expected_description


@pytest.mark.parametrize(
    "input_file, gene_data, do_report",
    [
        param("gene1_desc.json", {"gene_name": "gene1", "gene_desc": "old_desc"}, False, id="No report"),
        param("gene1_desc.json", {"gene_name": "gene1", "gene_desc": "old_desc"}, True, id="Do report"),
    ],
)
def test_load_description_do_report(
    annot_test_db: UnitTestDB,
    data_dir: Path,
    capsys: CaptureFixture,
    input_file: str,
    gene_data: dict[str, str],
    do_report: bool,
) -> None:
    """Tests the method `load_description()`"""
    with annot_test_db.dbc.test_session_scope() as session:
        add_gene(annot_test_db.dbc.dialect, session, gene_data)
        load_descriptions(session, data_dir / input_file, report=do_report)
        captured = capsys.readouterr()
        if do_report:
            assert str(captured.out).strip() == "gene\tgene1\tgene1\told_desc\tnew_desc"
        else:
            assert not captured.out
