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

import pytest
from pytest import param, raises

from ensembl.io.genomio.annotation import get_core_data, load_descriptions
from ensembl.core.models import Gene, Transcript, Translation, Xref, ObjectXref, metadata
from ensembl.utils.database import UnitTestDB


@pytest.fixture(name="test_db", scope="module")
def fixture_annotation_test_db(db_factory) -> UnitTestDB:
    """Returns a test database with all core tables and basic data."""
    test_db: UnitTestDB = db_factory(src="no_src", name="load_annotation")
    test_db.dbc.create_all_tables(metadata)
    # Add data: gene1 with 2 transcripts and translations, xref on gene1 only
    with test_db.dbc.session_scope() as session:
        session.add(
            Gene(
                gene_id=1,
                stable_id="gene1",
                biotype="protein_coding",
                description="gene1 description",
                analysis_id=1,
                seq_region_id=1,
                seq_region_start=1,
                seq_region_end=100,
                seq_region_strand=1,
                source="test_source",
                canonical_transcript_id=1,
            )
        )
        session.add(
            Transcript(
                gene_id=1,
                transcript_id=1,
                stable_id="gene1_t1",
                biotype="protein_coding",
                analysis_id=1,
                seq_region_id=1,
                seq_region_start=1,
                seq_region_end=100,
                seq_region_strand=1,
            )
        )
        session.add(
            Transcript(
                gene_id=1,
                transcript_id=2,
                stable_id="gene1_t2",
                biotype="protein_coding",
                analysis_id=1,
                seq_region_id=1,
                seq_region_start=1,
                seq_region_end=50,
                seq_region_strand=1,
            )
        )
        session.add(
            Translation(
                transcript_id=1,
                translation_id=1,
                stable_id="gene1_p1",
                seq_start=1,
                seq_end=100,
                start_exon_id=1,
                end_exon_id=2,
            )
        )
        session.add(
            Translation(
                transcript_id=2,
                translation_id=2,
                stable_id="gene1_p2",
                seq_start=1,
                seq_end=50,
                start_exon_id=1,
                end_exon_id=1,
            )
        )
        session.add(
            Xref(xref_id=1, dbprimary_acc="XREF_ACC1", display_label="xref1", external_db_id=1,)
        )
        session.add(
            ObjectXref(
                object_xref_id=1,
                xref_id=1,
                ensembl_object_type="gene",
                ensembl_id=1,
            )
        )
        session.commit()
    return test_db


@pytest.mark.parametrize(
    "table, expected_ids, expectation",
    [
        param("gene", ["gene1", "xref_acc1"], no_raise()),
        param("transcript", ["gene1_t1", "gene1_t2"], no_raise()),
        param("foo", [], raises(ValueError, match="Table foo is not supported")),
    ],
)
def test_get_core_data(
    test_db: UnitTestDB, table: str, expected_ids: list[str], expectation: ContextManager
) -> None:
    """Tests the method get_core_data()"""
    with test_db.dbc.test_session_scope() as session:
        with expectation:
            feats = get_core_data(session, table)
            assert list(feats.keys()) == expected_ids


@pytest.mark.parametrize(
    "input_file, report, do_update, match_xrefs, expectation",
    [
        param("group.json", False, False, False, no_raise()),
        param("group.json", False, True, False, no_raise()),
        param("group.json", False, False, True, no_raise()),
        param("gene1_xref.json", False, True, False, no_raise(), id="Match gene id"),
        param("gene1_alt_synonym.json", False, True, False, no_raise(), id="Match synonym"),
        param("gene1_alt_xref.json", False, True, True, no_raise(), id="Match xref"),
        param("gene1_alt_xref.json", False, True, False, no_raise(), id="Don't match xref"),
        param("gene1_alt_xref_nomatch.json", False, True, True, no_raise(), id="Match xref, no match"),
    ],
)
def test_load_descriptions(test_db: UnitTestDB, data_dir: Path, input_file: str, report: bool, do_update: bool, match_xrefs: bool, expectation: ContextManager):
    """Tests the method get_core_data()"""
    with test_db.dbc.test_session_scope() as session:
        with expectation:
            load_descriptions(session, data_dir / input_file, report=report, do_update=do_update, match_xrefs=match_xrefs)
