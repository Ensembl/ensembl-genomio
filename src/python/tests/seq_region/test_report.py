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
"""Unit testing of `ensembl.io.genomio.seq_region.report` module."""

from pathlib import Path

from pytest import raises

from ensembl.io.genomio.seq_region.report import ReportRecord


def test_from_report(data_dir: Path) -> None:
    """Test for `ReportRecord`."""
    report_file = "report.txt"
    record = ReportRecord(data_dir / report_file)
    seq = next(record.reader)

    expected_seq = {
        "Assembly-Unit": "Primary Assembly",
        "Assigned-Molecule": "Ia",
        "Assigned-Molecule-Location/Type": "Chromosome",
        "GenBank-Accn": "CM002034.1",
        "RefSeq-Accn": "NC_031467.1",
        "Relationship": "=",
        "Sequence-Length": "1859933",
        "Sequence-Name": "TGME49_chrIa",
        "Sequence-Role": "assembled-molecule",
        "UCSC-style-name": "na",
    }
    assert seq == expected_seq
    assert record.metadata.get("Assembly level") == "Chromosome"


def test_from_report_error(data_dir: Path) -> None:
    """Test for `ReportRecord` with a file without header."""
    report_file = "report_noheader.txt"
    with raises(ValueError):
        _ = ReportRecord(data_dir / report_file)
