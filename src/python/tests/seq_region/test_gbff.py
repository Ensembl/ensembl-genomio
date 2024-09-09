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
"""Unit testing of `ensembl.io.genomio.seq_region.gbff` module."""

from contextlib import nullcontext as no_raise
from pathlib import Path
from typing import ContextManager

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest
from pytest import param, raises

from ensembl.io.genomio.seq_region.gbff import GBFFRecord
from ensembl.io.genomio.seq_region.exceptions import UnknownMetadata


def get_record(gbff_path: Path) -> SeqRecord:
    """Returns the first SeqRecord from a Genbank file."""
    with gbff_path.open("r") as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            return record
    return SeqRecord(Seq(""))


def test_gbff(data_dir: Path) -> None:
    """Test for `GBFFRecord` init."""
    record = GBFFRecord(get_record(data_dir / "apicoplast.gb"))
    assert record


@pytest.mark.parametrize(
    "input_gb, expected_id",
    [
        param("apicoplast.gb", "U87145", id="Found genbank ID"),
        param("apicoplast_nocomment.gb", None, id="No comment"),
        param("apicoplast_simple_comment.gb", None, id="Comment without ID"),
    ],
)
def test_get_genbank_id(data_dir: Path, input_gb: str, expected_id: str | None) -> None:
    """Test for `gbff.get_genbank_id()`.

    Args:
        input_gb: Genbank file to use.
        expect_id: Expected Genbank ID.

    """
    seq_record = get_record(data_dir / input_gb)
    record = GBFFRecord(seq_record)
    assert record.get_genbank_id() == expected_id


@pytest.mark.parametrize(
    "input_gb, expected_table",
    [
        param("apicoplast.gb", 4, id="Found codon table"),
        param("apicoplast_nogene.gb", None, id="No codon table"),
    ],
)
def test_get_codon_table(data_dir: Path, input_gb: str, expected_table: str | None) -> None:
    """Test for `gbff.get_codon_table()`.

    Args:
        input_gb: Genbank file to use.
        expect_table: Expected codon table number.

    """
    seq_record = get_record(data_dir / input_gb)
    record = GBFFRecord(seq_record)
    assert record.get_codon_table() == expected_table


@pytest.mark.parametrize(
    "input_gb, expected_location, expectation",
    [
        param("apicoplast.gb", "apicoplast_chromosome", no_raise(), id="Found location"),
        param("apicoplast_nofeatures.gb", None, no_raise(), id="No features"),
        param("apicoplast_unknown_location.gb", None, raises(UnknownMetadata), id="Unknown location"),
        param("apicoplast_noprefix_location.gb", "apicoplast_chromosome", no_raise(), id="No prefix"),
    ],
)
def test_get_organelle(
    data_dir: Path, input_gb: str, expected_location: str | None, expectation: ContextManager
) -> None:
    """Test for `gbff.get_organelle()`.

    Args:
        input_gb: Genbank file to use.
        expect_location: Expected location of an organelle.

    """
    seq_record = get_record(data_dir / input_gb)
    record = GBFFRecord(seq_record)
    with expectation:
        assert record.get_organelle() == expected_location


def test_get_organelle_custom(data_dir: Path) -> None:
    """Test for `gbff.get_organelle()` with a custom location map."""
    input_gb = "apicoplast.gb"
    expected_location = "custom_location"
    custom_map = {"apicoplast": expected_location}
    seq_record = get_record(data_dir / input_gb)
    record = GBFFRecord(seq_record)
    assert record.get_organelle(custom_map) == expected_location


@pytest.mark.parametrize(
    "input_gb, expected_circular",
    [
        param("apicoplast.gb", True, id="Circular"),
        param("apicoplast_linear.gb", False, id="Not circular"),
    ],
)
def test_is_circular(data_dir: Path, input_gb: str, expected_circular: bool) -> None:
    """Test for `gbff.is_circular()`.

    Args:
        input_gb: Genbank file to use.
        expect_circular: If the sequence is circular.

    """
    seq_record = get_record(data_dir / input_gb)
    record = GBFFRecord(seq_record)
    assert record.is_circular() == expected_circular
