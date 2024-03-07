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
"""Unit testing of `ensembl.io.genomio.genome_metadata.prepare` module.

Typical usage example::
    $ pytest test_prepare.py

"""

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import Any, Callable, ContextManager, Dict, Optional
from unittest.mock import Mock, patch
from xml.etree import ElementTree

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_metadata import prepare


# @pytest.mark.dependency(name="test_get_gbff_regions")
@pytest.mark.parametrize(
    "genome_file, gff3_file, output, expectation",
    [
        pytest.param(
            "genbank_genome.json",
            None,
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
            "fake.gff3",
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
            "fake.gff3",
            {
                "assembly": {"provider_name": "GenBank", "provider_url": None},
                "annotation": {"provider_name": "GenBank", "provider_url": None},
            },
            does_not_raise(),
            id="Provider information already present",
        ),
        pytest.param(
            "cncb_genome.json", None, {}, pytest.raises(prepare.MetadataError), id="Unexpected provider"
        ),
    ],
)
def test_add_provider(
    json_data: Callable[[str], Any],
    genome_file: str,
    gff3_file: Optional[str],
    output: Dict[str, Dict[str, Optional[str]]],
    expectation: ContextManager,
) -> None:
    """Tests the `prepare.add_provider()` method.

    Args:
        json_data: JSON test file parsing fixture.
        genome_file: Genome metadata JSON file.
        gff3_file: GFF3 file.
        output: Expected elements present in the updated genome metadata.
        expectation: Context manager for the expected exception (if any).
    """
    genome_metadata = json_data(genome_file)
    with expectation:
        prepare.add_provider(genome_metadata, gff3_file)
        for section, metadata in output.items():
            for key, value in metadata.items():
                assert genome_metadata[section].get(key, None) == value


@pytest.mark.parametrize(
    "genome_file, output",
    [
        ("genbank_genome.json", 2),
        ("updated_genome.json", 1),
        ("cncb_genome.json", 0),
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
        ("genbank_genome.json", "03-2024"),
        ("updated_genome.json", "01-2021"),
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


@pytest.mark.dependency(name="test_get_node_text")
@pytest.mark.parametrize(
    "xml_file, tag, optional, output, expectation",
    [
        pytest.param(
            "",
            "tag",
            False,
            "",
            pytest.raises(prepare.MissingNodeError, match="No node provided to look for 'tag'"),
            id="No node provided",
        ),
        pytest.param("default_taxonomy.xml", "TAXON_ID", False, "34611", does_not_raise(), id="Tag present"),
        pytest.param("default_taxonomy.xml", "tag", True, None, does_not_raise(), id="Missing optional tag"),
        pytest.param(
            "default_taxonomy.xml",
            "tag",
            False,
            "",
            pytest.raises(prepare.MissingNodeError, match="No node found for tag 'tag'"),
            id="Missing mandatory tag",
        ),
    ],
)
def test_get_node_text(
    data_dir: Path,
    xml_file: str,
    tag: str,
    optional: bool,
    output: Optional[str],
    expectation: ContextManager,
) -> None:
    """Tests the `prepare._get_node_text()` method.

    Args:
        data_dir: Module's test data directory fixture.
        xml_file: XML file with assembly's taxonomy data.
        tag: Tag to fetch within the node.
        optional: Do not raise an exception if the tag does not exist.
        output: Expected field value returned.
        expectation: Context manager for the expected exception (if any).
    """
    if xml_file:
        tree = ElementTree.parse(data_dir / xml_file)
        node = tree.find(".//TAXON")
    else:
        node = None
    with expectation:
        result = prepare._get_node_text(node, tag, optional)
        assert result == output


@pytest.mark.dependency(depends=["test_get_node_text"])
@patch("requests.Response")
@patch("requests.get")
@pytest.mark.parametrize(
    "accession, base_api_url, xml_file, output, expectation",
    [
        pytest.param(
            "GCF_013436015.2",
            prepare.DEFAULT_API_URL,
            "default_taxonomy.xml",
            {"taxon_id": 34611, "scientific_name": "Rhipicephalus annulatus"},
            does_not_raise(),
            id="Basic taxonomy data",
        ),
        pytest.param(
            "GCA_013436015.2",
            "/",
            "strain_taxonomy.xml",
            {"taxon_id": 34611, "scientific_name": "Rhipicephalus annulatus", "strain": "Klein Grass"},
            does_not_raise(),
            id="Taxonomy with strain data",
        ),
        pytest.param(
            "GCA_013436015.2",
            "",
            "no_taxonomy.xml",
            {},
            pytest.raises(prepare.MissingNodeError, match="Cannot find the TAXON node"),
            id="Missing TAXON node",
        ),
    ],
)
def test_get_taxonomy_from_accession(
    mock_requests_get: Mock,
    mock_response: Mock,
    data_dir: Path,
    accession: str,
    base_api_url: str,
    xml_file: str,
    output: Dict[str, Any],
    expectation: ContextManager,
):
    """Tests the `prepare.get_taxonomy_from_accession()` method.

    Args:
        mock_requests_get: A mock of `requests.get()` function.
        mock_response: A mock of `requests.Response` class.
        data_dir: Module's test data directory fixture.
        accession: INSDC accession ID.
        base_api_url: Base API URL to fetch the taxonomy data from.
        xml_file: XML file with assembly's taxonomy data.
        output: Expected taxonomy data returned.
        expectation: Context manager for the expected exception (if any).
    """
    xml_path = data_dir / xml_file
    with xml_path.open() as xml:
        text = "".join(xml.readlines())
    mock_response.text = text
    mock_requests_get.return_value = mock_response
    with expectation:
        result = prepare.get_taxonomy_from_accession(accession, base_api_url)
        assert not DeepDiff(result, output)
