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
import pytest
from pytest import param

from ensembl.io.genomio.seq_region.gbff import GBFFRecord

def test_gbff(data_dir: Path):
    record = GBFFRecord(data_dir / "apicoplast.gb")
    assert record

def get_record(gbff_path: Path):
    with gbff_path.open("r") as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            return record

@pytest.mark.parametrize(
    "input_gb, expected_id",
    [
        param("apicoplast.gb", "U87145", id="Found genbank ID"),
    ],
)
def test_get_genbank_id(data_dir: Path, input_gb: str, expected_id: str | None):
    seq_record = get_record(data_dir / input_gb)
    record = GBFFRecord(seq_record)
    assert record.get_genbank_id() == expected_id
