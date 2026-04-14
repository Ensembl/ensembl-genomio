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
"""Unit testing of `ensembl.io.genomio.gff3.id_allocator` module."""

from contextlib import nullcontext as does_not_raise
from difflib import unified_diff
import filecmp
from pathlib import Path
from typing import ContextManager, Dict, List, Optional

import pytest

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from ensembl.io.genomio.gff3.id_allocator import StableIDAllocator, InvalidStableID
from ensembl.io.genomio.gff3.features import GFFSeqFeature
from pytest import raises


def _read_record(in_gff: Path) -> SeqRecord:
    """Get one clean record from a GFF3 file."""
    with in_gff.open("r") as in_gff_fh:
        for record in GFF.parse(in_gff_fh):
            new_record = SeqRecord(seq=record.seq, id=record.id)
            new_record.features = record.features
            return new_record

    raise ValueError(f"No feature loaded from {in_gff}")


def _write_record(out_record: SeqRecord, out_gff: Path) -> None:
    """Write one record to a GFF3 file."""
    with out_gff.open("w") as out_gff_fh:
        GFF.write([out_record], out_gff_fh)


def _show_diff(result_path: Path, expected_path: Path) -> str:
    """Create a useful diff between 2 files."""
    with open(result_path, "r") as result_fh:
        results = result_fh.readlines()
    with open(expected_path, "r") as expected_fh:
        expected = expected_fh.readlines()
    diff = list(unified_diff(expected, results))
    return "".join(diff)


@pytest.mark.parametrize(
    "genome, expected_prefix",
    [
        pytest.param({}, "TMP_PREFIX_", id="Default prefix"),
        pytest.param({"BRC4": {"organism_abbrev": "LOREM"}}, "TMP_LOREM_", id="Prefix from genome meta"),
    ],
)
def test_set_prefix(genome: Dict, expected_prefix: str) -> None:
    """Test prefix setting from genome metadata."""
    ids = StableIDAllocator()
    ids.set_prefix(genome)
    assert ids.prefix == expected_prefix


@pytest.mark.parametrize(
    "prefix, expected_ids",
    [
        pytest.param(None, ["TMP_1", "TMP_2"], id="Default prefix"),
        pytest.param("", ["1", "2"], id="No prefix as empty string"),
        pytest.param("MYPREF_", ["MYPREF_1", "MYPREF_2"], id="Prefix MYPREF_"),
    ],
)
def test_generate_id(prefix: str, expected_ids: List[str]) -> None:
    """Test IDs generation."""
    ids = StableIDAllocator()
    if prefix is not None:
        ids.prefix = prefix

    id1 = ids.generate_gene_id()
    id2 = ids.generate_gene_id()
    new_ids = [id1, id2]
    assert new_ids == expected_ids


@pytest.mark.parametrize(
    "min_id_length, test_id, outcome",
    [
        pytest.param(None, "LOREMIPSUM_01", True, id="OK ID"),
        pytest.param(None, "", False, id="Empty"),
        pytest.param(None, "A", False, id="Too short"),
        pytest.param(None, "Abc", False, id="Too short 2"),
        pytest.param(5, "Abcde", True, id="At custom min length"),
        pytest.param(5, "Abcd", False, id="Below custom min length"),
        pytest.param(None, "CHR1:100..200", False, id="Coordinates"),
        pytest.param(None, "LOREM|IPSUM", False, id="Special char |"),
        pytest.param(None, "LOREM IPSUM", False, id="Special char space"),
        pytest.param(None, "Trnaa-UAA", False, id="Trna ID, lower case"),
        pytest.param(None, "TRNAA-UAA", False, id="Trna ID, upper case"),
    ],
)
def test_valid_id(min_id_length: Optional[int], test_id: str, outcome: bool) -> None:
    """Test ID validity check."""
    ids = StableIDAllocator()
    if min_id_length is not None:
        ids.min_id_length = min_id_length

    assert ids.is_valid(test_id) == outcome


@pytest.mark.parametrize(
    "test_id, skip_flag, outcome",
    [
        pytest.param("LOREMIPSUM_01", True, True, id="Skip Good ID"),
        pytest.param("LO..rem|ipsum", True, True, id="Skip Bad ID"),
        pytest.param("LOREMIPSUM_01", False, True, id="No skip Good ID"),
        pytest.param("LO..rem|ipsum", False, False, id="No skip Bad ID"),
    ],
)
def test_valid_id_skip(test_id: str, skip_flag: bool, outcome: bool) -> None:
    """Test ID validity check without the validation flag."""
    ids = StableIDAllocator(skip_gene_id_validation=skip_flag)
    assert ids.is_valid(test_id) == outcome


@pytest.mark.parametrize(
    "test_id, prefixes, outcome",
    [
        pytest.param("LOREM-IPSUM1", [], "LOREM-IPSUM1", id="No prefixes"),
        pytest.param("LOREM-IPSUM1", ["DOLOR"], "LOREM-IPSUM1", id="Unused prefix"),
        pytest.param("LOREM-IPSUM1", ["LOREM-"], "IPSUM1", id="Found 1 prefix"),
        pytest.param("LOREM-IPSUM1", ["LOREM-", "IPSUM"], "IPSUM1", id="Only 1 prefix is removed"),
    ],
)
def test_remove_prefixes(test_id: str, prefixes: List[str], outcome: str) -> None:
    """Test prefix removal."""
    assert StableIDAllocator.remove_prefix(test_id, prefixes) == outcome


@pytest.mark.parametrize(
    "test_id, outcome",
    [
        pytest.param("LOREM-IPSUM1", "LOREM-IPSUM1", id="No prefixes"),
        pytest.param("cds-LOREM-IPSUM1", "LOREM-IPSUM1", id="Prefix cds-"),
        pytest.param("cds:LOREM-IPSUM1", "LOREM-IPSUM1", id="Prefix cds-"),
        pytest.param("bad", "", id="Short id without prefix"),
        pytest.param("cds:bad", "", id="Invalid id with cds:"),
    ],
)
def test_normalize_cds_id(test_id: str, outcome: str) -> None:
    """Test CDS id normalization."""
    ids = StableIDAllocator()
    assert ids.normalize_cds_id(test_id) == outcome


@pytest.mark.parametrize(
    "test_id, numbers, outcomes",
    [
        pytest.param("LOREM-IPSUM1", [1], ["LOREM-IPSUM1_t1"], id="1 transcript ID"),
        pytest.param("LOREM-IPSUM1", [1, 2], ["LOREM-IPSUM1_t1", "LOREM-IPSUM1_t2"], id="2 transcript IDs"),
        pytest.param("LOREM-IPSUM1", [1, 1], ["LOREM-IPSUM1_t1", "LOREM-IPSUM1_t1"], id="Same number (!)"),
    ],
)
def test_normalize_transcript_id(test_id: str, numbers: List[int], outcomes: List[str]) -> None:
    """Test transcript id normalization."""
    new_ids = []
    for number in numbers:
        new_ids.append(StableIDAllocator.generate_transcript_id(test_id, number))

    assert new_ids == outcomes


@pytest.mark.parametrize(
    "input_gff, expected_gff",
    [
        pytest.param("pseudo_01.gff3", "pseudo_01.gff3", id="Good ID, no change"),
        pytest.param("pseudo_02_in.gff3", "pseudo_02_out.gff3", id="invalid ID"),
        pytest.param("pseudo_03_in.gff3", "pseudo_03_out.gff3", id="Using gene ID"),
        pytest.param("pseudo_04.gff3", "pseudo_04.gff3", id="Without CDS"),
    ],
)
def test_normalize_pseudogene_cds_id(
    tmp_path: Path, data_dir: Path, input_gff: str, expected_gff: str
) -> None:
    """Test pseudogene CDS ID normalization."""
    ids = StableIDAllocator()

    # Load record and update feature
    record = _read_record(data_dir / input_gff)
    features = [GFFSeqFeature.cast(feat) for feat in record.features]
    record.features = []
    for feature in features:
        ids.normalize_pseudogene_cds_id(feature)
        record.features.append(feature)

    # Write it out and compare the GFF3 files
    out_path = tmp_path / "result.gff3"
    _write_record(record, out_path)
    expected_path = data_dir / expected_gff
    diff = _show_diff(out_path, expected_path)
    assert filecmp.cmp(out_path, expected_path), f"Files differ: {diff}"


@pytest.mark.parametrize(
    "input_gff, expected_id, make_id, expected",
    [
        pytest.param("geneid_ok.gff3", "LOREMIPSUM1", None, does_not_raise(), id="Good ID, no change"),
        pytest.param("geneid_makeid.gff3", "TMP_1", True, does_not_raise(), id="Make ID"),
        pytest.param("geneid_bad.gff3", "", False, raises(InvalidStableID), id="Bad ID, fail"),
        pytest.param("geneid_GeneID.gff3", "GeneID_000001", None, does_not_raise(), id="Replace with GeneID"),
        pytest.param("geneid_noGeneID.gff3", "TMP_1", True, does_not_raise(), id="Dbxref without Gene ID"),
    ],
)
def test_normalize_gene_id(
    data_dir: Path, input_gff: str, expected_id: str, make_id: Optional[bool], expected: ContextManager
) -> None:
    """Test gene ID normalization."""
    ids = StableIDAllocator()
    if make_id is not None:
        ids.make_missing_stable_ids = make_id

    # Load record and update feature
    record = _read_record(data_dir / input_gff)
    features = record.features
    record.features = []

    with expected:
        for feature in features:
            feature.id = ids.normalize_gene_id(GFFSeqFeature.cast(feature))
            assert feature.id == expected_id


@pytest.mark.parametrize(
    "input_gff, expected_ids, expected",
    [
        pytest.param(
            "geneid_GeneID2.gff3", ["GeneID_000001", "GeneID_000001_2"], does_not_raise(), id="Same GeneIDs"
        ),
    ],
)
def test_normalize_gene_id_duplicate(
    data_dir: Path, input_gff: str, expected_ids: List[str], expected: ContextManager
) -> None:
    """Test gene ID normalization with duplicate Gene ID features."""
    ids = StableIDAllocator()

    # Load record and update feature
    record = _read_record(data_dir / input_gff)
    features = record.features
    record.features = []

    found_ids = []
    with expected:
        for feature in features:
            new_feature_id = ids.normalize_gene_id(GFFSeqFeature.cast(feature))
            if new_feature_id and feature.id != new_feature_id:
                found_ids.append(new_feature_id)
        assert found_ids == expected_ids
