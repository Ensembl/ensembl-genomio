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
"""Unit testing of `ensembl.io.genomio.annotation.load` module.
"""

from contextlib import nullcontext as no_raise
from pathlib import Path
from typing import ContextManager
from sqlalchemy import text

import pytest
from pytest import param, raises

from ensembl.io.genomio.annotation import get_core_data, load_descriptions
from ensembl.core.models import Gene, Transcript, Translation, Xref, ObjectXref, metadata
from ensembl.utils.database import UnitTestDB


def add_gene(dialect: str, session, gene_data = dict) -> None:
    gene_name = gene_data.get("gene_name", "gene1")
    gene_description = gene_data.get("gene_description", "")
    tr_name = gene_data.get("tr_name", "tr1")
    tr_description = gene_data.get("tr_description", "")
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
    session.commit()
    if dialect == "mysql":
        session.execute(text("SET FOREIGN_KEY_CHECKS=1"))


@pytest.fixture(name="test_db", scope="module")
def fixture_annotation_test_db(db_factory) -> UnitTestDB:
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
    test_db: UnitTestDB, gene_data: dict, table: str, expected_ids: list[str], expectation: ContextManager
) -> None:
    """Tests the method get_core_data()"""
    with test_db.dbc.test_session_scope() as session:
        add_gene(test_db.dbc.dialect, session, gene_data)
        with expectation:
            feats = get_core_data(session, table)
            assert list(feats.keys()) == expected_ids


@pytest.mark.parametrize(
    "input_file, gene_data, table, expected_description",
    [
        param("gene1_nodesc.json", {"gene_name": "gene1"}, "gene", "", id="Gene: no desc -> no desc"),
        param("gene1_desc.json", {"gene_name": "gene1"}, "gene", "new_desc", id="Gene: no desc -> new desc"),
        param("gene1_desc.json", {"gene_name": "gene1", "gene_desc": "old_desc"}, "gene", "new_desc", id="Gene: old desc -> new desc"),
        param("gene1_desc.json", {"gene_name": "gene1", "gene_desc": "old_desc [Source: ext]"}, "gene", "new_desc", id="Gene: source desc -> new desc"),
        param("tr1_nodesc.json", {"tr_name": "tr1"}, "transcript", "", id="Tr: no desc -> no desc"),
        param("tr1_desc.json", {"tr_name": "tr1"}, "transcript", "new_desc", id="Tr: no desc -> new desc"),
        param("tr1_desc.json", {"tr_name": "tr1", "tr_desc": "old_desc"}, "transcript", "new_desc", id="Tr: old desc -> new desc"),
    ],
)
def test_load_description(
    test_db: UnitTestDB, data_dir: Path, input_file: str, gene_data: dict, table: str, expected_description: str) -> None:
    """Tests the method load_description()"""
    with test_db.dbc.test_session_scope() as session:
        add_gene(test_db.dbc.dialect, session, gene_data)
        load_descriptions(session, data_dir / input_file, do_update=True)
        feats = get_core_data(session, table)
        assert len(feats) == 1
        if table == "gene":
            name = gene_data["gene_name"]
        elif table == "transcript":
            name = gene_data["tr_name"]
        print(feats)
        assert feats[name]
        assert feats[name][2] == expected_description
