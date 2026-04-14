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
"""Unit testing of `ensembl.io.genomio.genome_metadata.prepare` module."""
# pylint: disable=too-many-positional-arguments

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import Any, Callable, ContextManager, Dict, Optional
from unittest.mock import Mock, patch

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_metadata import prepare


@pytest.mark.parametrize(
    "genome_file, ncbi_data, output, expectation",
    [
        pytest.param(
            "genbank_genome.json",
            {},
            {
                "assembly": {
                    "provider_name": "GenBank",
                    "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_013436015.2",
                },
            },
            does_not_raise(),
            id="GenBank assembly",
        ),
        pytest.param(
            "refseq_genome.json",
            {"annotation_info": {}},
            {
                "assembly": {
                    "provider_name": "RefSeq",
                    "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004695.1",
                },
                "annotation": {
                    "provider_name": "RefSeq",
                    "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004695.1",
                },
            },
            does_not_raise(),
            id="RefSeq assembly and annotation",
        ),
        pytest.param(
            "updated_genome.json",
            {"annotation_info": {}},
            {
                "assembly": {"provider_name": "GenBank", "provider_url": None},
                "annotation": {"provider_name": "GenBank", "provider_url": None},
            },
            does_not_raise(),
            id="Provider information already present",
        ),
        pytest.param(
            "cncb_genome.json", {}, {}, pytest.raises(prepare.MetadataError), id="Unexpected provider"
        ),
    ],
)
def test_add_provider(
    json_data: Callable[[str], Any],
    genome_file: str,
    ncbi_data: Dict,
    output: Dict[str, Dict[str, Optional[str]]],
    expectation: ContextManager,
) -> None:
    """Tests the `prepare.add_provider()` method.

    Args:
        json_data: JSON test file parsing fixture.
        genome_file: Genome metadata JSON file.
        ncbi_data: Report from NCBI datasets.
        output: Expected elements present in the updated genome metadata.
        expectation: Context manager for the expected exception (if any).
    """
    genome_metadata = json_data(genome_file)
    with expectation:
        prepare.add_provider(genome_metadata, ncbi_data)
        for section, metadata in output.items():
            for key, value in metadata.items():
                assert genome_metadata[section].get(key, None) == value


@pytest.mark.parametrize(
    "genome_file, output",
    [
        pytest.param("genbank_genome.json", 2, id="Added assembly version"),
        pytest.param("updated_genome.json", 1, id="Version found, nothing to add"),
        pytest.param("cncb_genome.json", 0, id="No version available, nothing to add"),
    ],
)
def test_add_assembly_version(json_data: Callable[[str], Any], genome_file: str, output: int) -> None:
    """Tests the `prepare.add_assembly_version()` method.

    Args:
        json_data: JSON test file parsing fixture.
        genome_file: Genome metadata JSON file.
        output: Assembly version expected in the updated genome metadata.
    """
    genome_metadata = json_data(genome_file)
    prepare.add_assembly_version(genome_metadata)
    assert genome_metadata["assembly"].get("version", 0) == output


@patch("datetime.date")
@pytest.mark.parametrize(
    "genome_file, output",
    [
        pytest.param("genbank_genome.json", "03-2024", id="Added '03-2024' as genebuild metadata"),
        pytest.param("updated_genome.json", "01-2021", id="Found '01-2021', nothing to add"),
    ],
)
def test_add_genebuild_metadata(
    mock_date: Mock, json_data: Callable[[str], Any], genome_file: str, output: str
) -> None:
    """Tests the `prepare.add_genebuild_metadata()` method.

    Args:
        mock_date: A mock of `datetime.date` class.
        json_data: JSON test file parsing fixture.
        genome_file: Genome metadata JSON file.
        output: Expected date for genebuild's `"start_date"` and `"version"` in the updated genome metadata.
    """
    mock_date.today.return_value = mock_date
    mock_date.isoformat.return_value = output
    genome_metadata = json_data(genome_file)
    prepare.add_genebuild_metadata(genome_metadata)
    assert genome_metadata["genebuild"]["start_date"] == output
    assert genome_metadata["genebuild"]["version"] == output


@pytest.mark.parametrize(
    "genome_file, ncbi_data_organism, output",
    [
        pytest.param(
            "genbank_genome.json",
            {"tax_id": 34611, "organism_name": "Rhipicephalus annulatus"},
            {"taxonomy_id": 34611, "scientific_name": "Rhipicephalus annulatus"},
            id="Add taxonomy information",
        ),
        pytest.param(
            "refseq_genome.json",
            {
                "tax_id": 34611,
                "organism_name": "Rhipicephalus annulatus",
                "infraspecific_names": {"strain": "Klein Grass"},
            },
            {"taxonomy_id": 34611, "scientific_name": "Rhipicephalus annulatus", "strain": "Klein Grass"},
            id="Add strain taxonomy information",
        ),
        pytest.param(
            "updated_genome.json",
            {"tax_id": 34611},
            {"taxonomy_id": 10092, "scientific_name": "Mus musculus", "strain": "domesticus"},
            id="Nothing to add",
        ),
    ],
)
def test_add_species_metadata(
    json_data: Callable[[str], Any],
    genome_file: str,
    ncbi_data_organism: Dict,
    output: Dict[str, Any],
) -> None:
    """Tests the `prepare.add_species_metadata()` method.

    Args:
        json_data: JSON test file parsing fixture.
        genome_file: Genome metadata JSON file.
        ncbi_data_organism: NCBI dataset organism report.
        output: Expected `"species"` genome metadata content.
    """
    ncbi_data = {"organism": ncbi_data_organism}
    genome_metadata = json_data(genome_file)
    prepare.add_species_metadata(genome_metadata, ncbi_data)
    assert not DeepDiff(genome_metadata["species"], output)


@patch("datetime.date")
@pytest.mark.parametrize(
    "input_filename, ncbi_filename, expected_filename",
    [
        pytest.param(
            "genome_accession.json",
            "genome_accession_ncbi.json",
            "genome_accession_updated.json",
            id="Add information from an accession",
        ),
        pytest.param(
            "updated_genome.json",
            "updated_genome_ncbi.json",
            "updated_genome.json",
            id="Nothing to add",
        ),
    ],
)
def test_prepare_genome_metadata(
    mock_date: Mock,
    tmp_path: Path,
    data_dir: Path,
    assert_files: Callable[[Path, Path], None],
    input_filename: str,
    ncbi_filename: str,
    expected_filename: str,
) -> None:
    """Tests the `prepare.prepare_genome_metadata()` method.

    Args:
        input_filename: Input genome JSON file.
        ncbi_filename: NCBI dataset report JSON file.
        expected_filename: Expected output genome JSON file.
    """
    mock_date.today.return_value = mock_date
    mock_date.isoformat.return_value = "2024-03-19"

    input_file = data_dir / input_filename
    ncbi_meta = data_dir / ncbi_filename
    output_file = tmp_path / expected_filename
    expected_file = data_dir / expected_filename
    prepare.prepare_genome_metadata(input_file=input_file, output_file=output_file, ncbi_meta=ncbi_meta)
    assert_files(output_file, expected_file)
