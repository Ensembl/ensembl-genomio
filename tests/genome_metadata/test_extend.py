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
"""Unit testing of `ensembl.io.genomio.genome_metadata.extend` module."""
# pylint: disable=too-many-positional-arguments

from pathlib import Path
from typing import Callable, Dict, List, Tuple

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_metadata import extend


@pytest.mark.dependency(name="test_get_gbff_regions")
@pytest.mark.parametrize(
    "gbff_file, output",
    [
        pytest.param("", [], id="No GBFF file"),
        pytest.param("sequences.gbff", ["CP089274", "CP089275", "RefChr0002"], id="sequences.gbff"),
    ],
)
def test_get_gbff_regions(data_dir: Path, gbff_file: str, output: List[str]) -> None:
    """Tests the `extend.get_gbff_regions()` method.

    Args:
        data_dir: Module's test data directory fixture.
        gbff_file: GBFF file name.
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
                {"Assembly name": "GCA000002765", "Organism name": "Plasmodium falciparum", "Taxid": "36329"},
            ),
            id="assembly_report.txt",
        ),
    ],
)
def test_report_to_csv(data_dir: Path, report_file: str, output: Tuple[str, Dict]) -> None:
    """Tests the `extend._report_to_csv()` method.

    Args:
        data_dir: Module's test data directory fixture.
        report_file: Assembly report file name.
        output: Expected returned value for the given assembly report file.
    """
    report_path = data_dir / report_file
    # pylint: disable=protected-access
    result = extend._report_to_csv(report_path)
    assert result[0] == output[0]
    assert not DeepDiff(result[1], output[1])


@pytest.mark.dependency(name="test_get_report_regions_names", depends=["test_report_to_csv"])
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
    """Tests the `extend.get_report_regions_names()` method.

    Args:
        data_dir: Module's test data directory fixture.
        report_file: Assembly report file name.
        output: Expected returned value for the given assembly report file.
    """
    report_path = data_dir / report_file
    result = extend.get_report_regions_names(report_path)
    assert result == output


@pytest.mark.dependency(
    name="test_get_additions", depends=["test_get_gbff_regions", "test_get_report_regions_names"]
)
@pytest.mark.parametrize(
    "report_file, gbff_file, output",
    [
        pytest.param(
            "assembly_report.txt", "", ["CP089275", "RefChr0001", "RefChr0002"], id="Additional regions found"
        ),
        pytest.param("assembly_report.txt", "sequences.gbff", [], id="No additional regions"),
    ],
)
def test_get_additions(data_dir: Path, report_file: str, gbff_file: str, output: List[str]) -> None:
    """Tests the `extend.get_additions()` method.

    Args:
        data_dir: Module's test data directory fixture.
        report_file: Assembly report file name.
        gbff_path: GBFF file name.
        output: Expected sequence regions names that need to be added.
    """
    report_path = data_dir / report_file
    gbff_path = data_dir / gbff_file if gbff_file else None
    result = extend.get_additions(report_path, gbff_path)
    assert result == output


@pytest.mark.dependency(depends=["test_get_additions"])
@pytest.mark.parametrize(
    "genome_infile, report_file, genbank_file, output_file",
    [
        pytest.param("genome.json", "", "", "genome.json", id="No report file"),
        pytest.param(
            "genome.json", "assembly_report.txt", "", "updated_genome.json", id="Additional seq regions"
        ),
        pytest.param(
            "genome.json", "assembly_report.txt", "sequences.gbff", "genome.json", id="No additional regions"
        ),
    ],
)
def test_amend_genome_metadata(
    tmp_path: Path,
    data_dir: Path,
    assert_files: Callable[[Path, Path], None],
    genome_infile: str,
    report_file: str,
    genbank_file: str,
    output_file: str,
) -> None:
    """Tests the `extend.amend_genome_metadata()` method.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        data_dir: Module's test data directory fixture.
        assert_files: File diff assertion fixture.
        genome_infile: Input genome metadata file.
        report_file: INSDC/RefSeq sequences report file.
        genbank_file: INSDC/RefSeq GBFF file.
        output_file: Expected amended genome metadata file.
    """
    genome_inpath = data_dir / genome_infile
    report_path = data_dir / report_file if report_file else None
    genbank_path = data_dir / genbank_file if genbank_file else None
    genome_outpath = tmp_path / "genome.out"
    extend.amend_genome_metadata(genome_inpath, genome_outpath, report_path, genbank_path)
    assert_files(genome_outpath, data_dir / output_file)
