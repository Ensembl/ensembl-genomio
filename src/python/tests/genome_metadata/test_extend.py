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
"""Unit testing of `ensembl.io.genomio.genome_metadata.extend` module.

Typical usage example::
    $ pytest test_extend.py

"""

from pathlib import Path
from typing import Dict, List, Tuple

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_metadata import extend


@pytest.mark.parametrize(
    "gbff_file, output",
    [
        pytest.param("", [], id="No GBFF file"),
        ("input.gbff.gz", ["LR605957", "LR605956"]),
    ],
)
def test_get_gbff_regions(data_dir: Path, gbff_file: str, output: List[str]) -> None:
    """Tests the `extend.get_gbff_regions` class.

    Args:
        data_dir: Module's test data directory fixture.
        gbff_path: GBFF file name.
        output: Expected list of sequence region IDs.
    """
    if gbff_file:
        gbff_path = data_dir / gbff_file
    else:
        gbff_path = None
    result = extend.get_gbff_regions(gbff_path)
    assert result == output


@pytest.mark.dependency(name="test_report_to_csv")
@pytest.mark.parametrize(
    "report_file, output",
    [
        pytest.param(
            "no_metadata_report.txt",
            ("1\t1\tChromosome\tCP089274.1\tRefChr0001.1\t5935961", {}),
            id="no_metadata_report.txt",
        ),
        pytest.param(
            "assembly_report.txt",
            (
                (
                    "Name\tMolecule\tLocation\tGenBank-Accn\tRefSeq-Accn\tLength\n"
                    "1\t1\tChromosome\tCP089274.1\tRefChr0001.1\t5935961\n"
                    "2\t2\tChromosome\tCP089275.1\tna\t5880203\n3\t3\tChromosome\tna\tRefChr0002.1\t5901247"
                ),
                {"Assembly name": "ASM2392016v1", "Organism name": "Curvularia clavata", "Taxid": "95742"},
            ),
            id="assembly_report.txt",
        ),
    ],
)
def test_report_to_csv(data_dir: Path, report_file: str, output: Tuple[str, Dict]) -> None:
    """Tests the `extend._report_to_csv` class.

    Args:
        data_dir: Module's test data directory fixture.
        report_file: Assembly report file name.
        output: Expected returned value for the given assembly report file.
    """
    report_path = data_dir / report_file
    result = extend._report_to_csv(report_path)
    assert result[0] == output[0]
    assert not DeepDiff(result[1], output[1])


@pytest.mark.dependency(depends=["test_report_to_csv"])
@pytest.mark.parametrize(
    "report_file, output",
    [
        pytest.param(
            "assembly_report.txt",
            [("CP089274", "RefChr0001"), ("CP089275", ""), ("", "RefChr0002")],
            id="assembly_report.txt",
        ),
    ],
)
def test_get_report_regions_names(data_dir: Path, report_file: str, output: List[Tuple[str, str]]) -> None:
    """Tests the `extend.get_report_regions_names` class.

    Args:
        data_dir: Module's test data directory fixture.
        report_file: Assembly report file name.
        output: Expected returned value for the given assembly report file.
    """
    report_path = data_dir / report_file
    result = extend.get_report_regions_names(report_path)
    assert result == output
