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
"""Unit testing of `ensembl.io.genomio.gff3.simplifier` Record object."""

from contextlib import nullcontext as does_not_raise
from os import PathLike
from pathlib import Path
from typing import Any, Callable, ContextManager, Dict, List, Optional, Union

from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.simplifier import Records


@pytest.mark.parametrize(
        "in_gff, excluded, expected_loaded",
        [
            param("record_n1.gff", None, ["scaffold1"], id="1 record"),
            param("record_n2.gff", None, ["scaffold1", "scaffold2"], id="2 records"),
            param("record_n2.gff", ["scaffold1"], ["scaffold2"], id="2 records, exclude 1"),
            param("record_n1.gff", ["Lorem"], ["scaffold1"], id="1 record, exclude not in record"),
        ],
)
def test_from_gff(tmp_dir: Path, data_dir: Path, in_gff: PathLike, excluded: Optional[List[str]], expected_loaded: List[str]) -> None:
    input = data_dir / in_gff

    records = Records()
    if excluded is None:
        records.from_gff(input)
    else:
        records.from_gff(input, excluded)
    record_names = [record.id for record in records]
    assert record_names == expected_loaded


@pytest.mark.parametrize(
        "in_gff",
        [
            param("record_n1.gff", id="1 record"),
            param("record_n2.gff", id="2 records"),
        ],
)
def test_to_gff(tmp_dir: Path, data_dir: Path, assert_files: Callable, in_gff: PathLike) -> None:
    input = data_dir / in_gff
    output = tmp_dir / in_gff
    records = Records()
    records.from_gff(input)
    records.to_gff(output)
    assert_files(input, output)
