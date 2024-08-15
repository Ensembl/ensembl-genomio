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
from typing import Any, ContextManager
from unittest.mock import MagicMock

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pytest
from pytest import param, raises

from ensembl.io.genomio.seq_region.exceptions import UnknownMetadata
from ensembl.io.genomio.seq_region.gbff import GBFFRecord
from ensembl.io.genomio.seq_region.collection import SeqCollection


def get_record(gbff_path: Path) -> SeqRecord:
    """Returns the first SeqRecord from a Genbank file."""
    with gbff_path.open("r") as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            return record


def test_collection():
    """Test for `SeqCollection` init."""
    seqs = SeqCollection()
    assert not seqs.seqs


@pytest.mark.parametrize(
    "genbank_id, codon_table, location, is_circular, expected_seq",
    [
        param(None, None, None, False, {"length": 3}, id="No extra data"),
        param(
            "Foo",
            None,
            None,
            False,
            {"length": 3, "synonyms": [{"source": "INSDC", "name": "Foo"}]},
            id="With genbank ID",
        ),
        param(None, None, None, True, {"length": 3, "circular": True}, id="With circular"),
        param(None, 4, None, False, {"length": 3, "codon_table": 4}, id="With codon table"),
        param(None, None, "apicoplast", False, {"length": 3, "location": "apicoplast"}, id="With location"),
    ],
)
def test_make_seqregion_from_gbff(
    genbank_id: str | None,
    codon_table: int | None,
    location: str | None,
    is_circular: bool,
    expected_seq: dict,
):
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


def test_from_gbff(data_dir: Path):
    """Test for `SeqCollection.from_gbff`."""
    gb_file = "apicoplast.gb"
    expected_seqs = {
        "NC_001799.1": {
            "circular": True,
            "codon_table": 4,
            "length": 60,
            "location": "apicoplast_chromosome",
            "synonyms": [{"name": "U87145", "source": "INSDC"}],
        }
    }
    collection = SeqCollection()
    collection.from_gbff(data_dir / gb_file)
    assert collection.seqs == expected_seqs


base_seq = {
    "Assembly-Unit": "Primary Assembly",
    "Assigned-Molecule": "I",
    "Assigned-Molecule-Location/Type": "Chromosome",
    "GenBank-Accn": "gb_id",
    "RefSeq-Accn": "rs_id",
    "Relationship": "=",
    "Sequence-Length": "1000",
    "Sequence-Name": "seq_name",
    "Sequence-Role": "assembled-molecule",
    "UCSC-style-name": "na",
}


@pytest.mark.parametrize(
    "seq_data, is_refseq, expected_key, expected_value, expected",
    [
        param({"GenBank-Accn": "na"}, False, None, None, no_raise(), id="Missing Genbank ID"),
        param({"RefSeq-Accn": "na"}, True, None, None, no_raise(), id="Missing RefSeq ID"),
        param({}, False, "name", "gb_id", no_raise(), id="use Genbank ID"),
        param({}, True, "name", "rs_id", no_raise(), id="use RefSeq ID"),
        param(
            {"Assigned-Molecule": "na", "GenBank-Accn": "na", "Sequence-Name": "na"},
            False,
            "synonyms",
            None,
            no_raise(),
            id="No synonyms",
        ),
        param(
            {"Assigned-Molecule": "na", "RefSeq-Accn": "na", "GenBank-Accn": "Foo", "Sequence-Name": "na"},
            False,
            "synonyms",
            [{"source": "GenBank", "name": "Foo"}],
            no_raise(),
            id="1 synonym",
        ),
        param({"Sequence-Length": "na"}, False, "length", None, no_raise(), id="No length"),
        param(
            {"Sequence-Role": "unplaced-scaffold"},
            False,
            "coord_system_level",
            "scaffold",
            no_raise(),
            id="Unplaced_scaffold level",
        ),
        param(
            {"Sequence-Role": "assembled-molecule"},
            False,
            "coord_system_level",
            "chromosome",
            no_raise(),
            id="Chromosome level",
        ),
        param(
            {"Sequence-Role": "assembled-molecule", "Assigned-Molecule-Location/Type": "plasmid"},
            False,
            "location",
            "plasmid",
            no_raise(),
            id="Chromosome level, plasmid location",
        ),
        param({"Sequence-Role": "foo"}, False, None, None, raises(UnknownMetadata), id="Unsupported role"),
        param(
            {"Sequence-Role": "assembled-molecule", "Assigned-Molecule-Location/Type": "foo"},
            False,
            None,
            None,
            raises(UnknownMetadata),
            id="Unsupported location",
        ),
    ],
)
def test_make_seqregion_from_report(
    seq_data: dict, is_refseq: bool, expected_key: str, expected_value: Any | None, expected: ContextManager
):
    collection = SeqCollection()
    input_data = {}
    input_data.update(base_seq)
    input_data.update(seq_data)
    with expected:
        seq_dict = collection.make_seq_region_from_report(input_data, is_refseq=is_refseq)
        if expected_key:
            assert seq_dict.get(expected_key) == expected_value
        else:
            assert not seq_dict

def test_remove():
    """Test `SeqCollection.remove()`."""
    collection = SeqCollection()
    collection.seqs = {"Foo": {}, "Bar": {}, "Lorem": {}}
    collection.remove(["Foo", "Lorem"])
    assert collection.seqs == {"Bar": {}}
    collection.remove(["OTHER"])
    assert collection.seqs == {"Bar": {}}


@pytest.mark.parametrize(
    "input_seq, code_map, expected_codon_table",
    [
        param({}, None, None, id="No data, no codon table"),
        param({"codon_table": 4}, None, 4, id="Existing codon_table"),
        param({"location": "foo"}, {"foo": 4}, 4, id="Location, new codon table"),
        param({"location": "foo"}, {"bar": 4}, None, id="Unsupported location, no table"),
        param({"codon_table": 2, "location": "foo"}, {"foo": 4}, 2, id="Existing codon table, keep"),
        param({"location": "apicoplast_chromosome"}, None, 4, id="use default map"),
    ],
)
def test_add_translation_table(input_seq: dict[str, str], code_map: dict[str, int] | None, expected_codon_table: int | None):
    collection = SeqCollection()
    seqname = "seq_name"
    collection.seqs = {seqname: input_seq}
    collection.add_translation_table(code_map)
    assert collection.seqs["seq_name"].get("codon_table") == expected_codon_table
