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

from typing import Optional

import pytest

from ensembl.io.genomio.annotation import get_core_data, load_descriptions
from ensembl.core.models import Gene, Transcript, Translation, metadata
from ensembl.utils.database import UnitTestDB


@pytest.fixture(name="test_db", scope="module")
def fixture_annotation_test_db(db_factory) -> UnitTestDB:
    """Returns a test database with all core tables and basic data."""
    test_db: UnitTestDB = db_factory("", "load_annotation")
    test_db.dbc.create_all_tables(metadata)
    # Add data
    with test_db.dbc.session_scope() as session:
        session.add(
            Gene(
                gene_id=1,
                stable_id="gene1",
                biotype="protein_coding",
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
        session.commit()
    return test_db


def test_get_core_data(test_db: UnitTestDB) -> None:
    """Tests the method get_core_data()"""
    with test_db.dbc.test_session_scope() as session:
        genes = get_core_data(session, "gene")
        assert list(genes.keys()) == ["gene1"]
        transcripts = get_core_data(session, "transcript")
        assert list(transcripts.keys()) == ["gene1_t1", "gene1_t2"]
