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
"""Unit testing of `ensembl.io.genomio.seq_region.collection` module."""

from contextlib import nullcontext as no_raise
from pathlib import Path
from typing import ContextManager

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pytest
from pytest import param, raises
from unittest.mock import MagicMock

from ensembl.io.genomio.seq_region.gbff import GBFFRecord
from ensembl.io.genomio.seq_region.collection import SeqCollection
from ensembl.io.genomio.seq_region.exceptions import UnknownMetadata


def get_record(gbff_path: Path) -> SeqRecord:
    """Returns the first SeqRecord from a Genbank file."""
    with gbff_path.open("r") as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            return record


def test_collection():
    """Test for `SeqCollection` init."""
    seqs = SeqCollection()
    assert seqs.seqs == {}


@pytest.mark.parametrize(
    "genbank_id, codon_table, location, is_circular, expected_seq",
    [
        param(None, None, None, False, {"length": 3}, id="No extra data"),
        param("Foo", None, None, False, {"length": 3, "synonyms": [{"source": "INSDC", "name": "Foo"}]}, id="With genbank ID"),
        param(None, None, None, True, {"length": 3, "circular": True}, id="With circular"),
        param(None, 4, None, False, {"length": 3, "codon_table": 4}, id="With codon table"),
        param(None, None, "apicoplast", False, {"length": 3, "location": "apicoplast"}, id="With location"),
    ],
)
def test_make_seqregion_from_gbff(genbank_id: str | None, codon_table: int | None, location: str | None, is_circular: bool, expected_seq: dict):
    """Test for `SeqCollection.make_seq_dict()`.

    Args:
        input_gb: Genbank file to use.
        expect_id: Expected Genbank ID.

    """
    record = SeqRecord(id="test_seq", seq=Seq("AAA"))
    record_data = GBFFRecord(record)
    record_data.get_genbank_id = MagicMock(return_value=genbank_id)
    record_data.get_codon_table = MagicMock(return_value=codon_table)
    record_data.get_organelle = MagicMock(return_value=location)
    record_data.is_circular = MagicMock(return_value=is_circular)
    collection = SeqCollection()
    seq_dict = collection.make_seqregion_from_gbff(record_data)
    assert seq_dict == expected_seq
