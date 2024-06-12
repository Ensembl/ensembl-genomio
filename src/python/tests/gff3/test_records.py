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
"""Unit testing of `ensembl.io.genomio.gff3.simplifier.Records` class."""

from contextlib import nullcontext as no_raise
from os import PathLike
from pathlib import Path
from typing import Callable, ContextManager, List, Optional

import pytest
from pytest import param, raises

from ensembl.io.genomio.gff3.simplifier import Records


@pytest.mark.parametrize(
    "in_gff, excluded, expected_loaded, expectation",
    [
        param("record_n2.gff", None, ["scaffold1", "scaffold2"], no_raise(), id="2 records"),
        param("record_n2.gff", ["scaffold1"], ["scaffold2"], no_raise(), id="2 records, exclude 1"),
        param("record_n1.gff", ["Lorem"], ["scaffold1"], no_raise(), id="1 record, exclude not in record"),
        param("invalid.gff", None, [], raises(AssertionError), id="Invalid GFF3"),
    ],
)
def test_from_gff(
    data_dir: Path,
    in_gff: PathLike,
    excluded: Optional[List[str]],
    expected_loaded: List[str],
    expectation: ContextManager,
) -> None:
    """Test loading GFF records from file."""
    input_gff = data_dir / in_gff

    records = Records()
    with expectation:
        records.from_gff(input_gff, excluded)
    if expected_loaded:
        record_names = [record.id for record in records]
        assert record_names == expected_loaded


@pytest.mark.parametrize(
    "in_gff",
    [
        param("record_n1.gff", id="1 record"),
        param("record_n2.gff", id="2 records"),
    ],
)
def test_to_gff(tmp_path: Path, data_dir: Path, assert_files: Callable, in_gff: PathLike) -> None:
    """Test writing GFF records to file."""
    input_gff = data_dir / in_gff
    output_gff = tmp_path / in_gff
    records = Records()
    records.from_gff(input_gff)
    records.to_gff(output_gff)
    assert_files(input_gff, output_gff)
