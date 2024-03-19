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
"""Expand the genome metadata file adding information about the provider, taxonomy, and assembly and
gene build versions.
"""

__all__ = [
    "add_provider",
    "add_assembly_version",
    "add_genebuild_metadata",
    "add_species_metadata",
    "get_taxonomy_from_accession",
    "prepare_genome_metadata",
    "PROVIDER_DATA",
    "DEFAULT_API_URL",
    "MissingNodeError",
    "MetadataError",
]

import datetime
from os import PathLike
from typing import Any, Dict, Optional
from xml.etree import ElementTree
from xml.etree.ElementTree import Element

import requests

from ensembl.io.genomio.utils import get_json, print_json
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


PROVIDER_DATA = {
    "GenBank": {
        "assembly": {
            "provider_name": "GenBank",
            "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome",
        },
        "annotation": {
            "provider_name": "GenBank",
            "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome",
        },
    },
    "RefSeq": {
        "assembly": {
            "provider_name": "RefSeq",
            "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome",
        },
        "annotation": {
            "provider_name": "RefSeq",
            "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome",
        },
    },
}
DEFAULT_API_URL = "https://www.ebi.ac.uk/ena/browser/api/xml"


class MissingNodeError(Exception):
    """When a taxon XML node cannot be found."""


class MetadataError(Exception):
    """When a metadata value is not expected."""


def add_provider(genome_metadata: Dict, gff3_file: Optional[PathLike] = None) -> None:
    """Updates the genome metadata adding provider information for assembly and gene models.

    Assembly provider metadata will only be added if it is missing, i.e. neither `"provider_name"` or
    `"provider_url"` are present. The gene model metadata will only be added if `gff3_file` is provided.

    Args:
        genome_data: Genome information of assembly, accession and annotation.
        gff3_file: Path to GFF3 file to use as annotation source for this genome.

    Raises:
        MetadataError: If accession's format in genome metadata does not match with a known provider.
    """
    # Get accession provider
    accession = genome_metadata["assembly"]["accession"]
    if accession.startswith("GCF"):
        provider = PROVIDER_DATA["RefSeq"]
    elif accession.startswith("GCA"):
        provider = PROVIDER_DATA["GenBank"]
    else:
        raise MetadataError(f"Accession does not look like an INSDC or RefSeq accession: {accession}")

    # Add assembly provider (if missing)
    assembly = genome_metadata["assembly"]
    if (not "provider_name" in assembly) and (not "provider_url" in assembly):
        assembly["provider_name"] = provider["assembly"]["provider_name"]
        assembly["provider_url"] = f'{provider["assembly"]["provider_url"]}/{accession}'

    # Add annotation provider if there are gene models
    if gff3_file:
        annotation = genome_metadata.setdefault("annotation", {})
        if ("provider_name" not in annotation) and ("provider_url" not in annotation):
            annotation["provider_name"] = provider["annotation"]["provider_name"]
            annotation["provider_url"] = f'{provider["annotation"]["provider_url"]}/{accession}'


def add_assembly_version(genome_data: Dict) -> None:
    """Adds version number to the genome's assembly information if one is not present already.

    Args:
        genome_data: Genome information of assembly, accession and annotation.
    """
    assembly = genome_data["assembly"]
    if not "version" in assembly:
        accession = assembly["accession"]
        version = accession.partition(".")[2]
        if version:
            assembly["version"] = int(version)


def add_genebuild_metadata(genome_data: Dict) -> None:
    """Adds genebuild metadata to genome information if not present already.

    The default convention is to use the current date as `"version"` and `"start_date"`.

    Args:
        genome_data: Genome information of assembly, accession and annotation.
    """
    genebuild = genome_data.setdefault("genebuild", {})
    current_date = datetime.date.today().isoformat()
    if not "version" in genebuild:
        genebuild["version"] = current_date
    if not "start_date" in genebuild:
        genebuild["start_date"] = current_date


def add_species_metadata(genome_metadata: Dict, base_api_url: str = DEFAULT_API_URL) -> None:
    """Adds missing species metadata from its taxonomy based on the genome's accession.

    If `"taxonomy_id"` is already present in the species metadata, nothing is added. The `"taxonomy_id"`,
    `"scientific_name"` and `"strain"` will be fetched from the taxonomy information linked to the given
    accession.

    Args:
        genome_metadata: Genome information of assembly, accession and annotation.
        base_api_url: Base API URL to fetch the taxonomy data from.

    """
    species = genome_metadata.setdefault("species", {})
    if not "taxonomy_id" in species:
        accession = genome_metadata["assembly"]["accession"]
        taxonomy = get_taxonomy_from_accession(accession, base_api_url)
        species["taxonomy_id"] = taxonomy["taxon_id"]
        if (not "strain" in species) and ("strain" in taxonomy):
            species["strain"] = taxonomy["strain"]
        if not "scientific_name" in species:
            species["scientific_name"] = taxonomy["scientific_name"]


def get_taxonomy_from_accession(accession: str, base_api_url: str = DEFAULT_API_URL) -> Dict[str, Any]:
    """Returns the taxonomy metadata associated to the given accession.

    Args:
        accession: INSDC accession ID.
        base_api_url: Base API URL to fetch the taxonomy data from.

    Returns:
        Dictionary with key-value pairs for `"taxon_id"` and `"scientific_name"`. `"strain"` will also be
        included if it is present in the fetched taxonomy data.

    Raises:
        MissingNodeError: If `"TAXON"` node is missing in the taxonomy data fetched.

    """
    # Use the GenBank accession without version
    gb_accession = accession.replace("GCF", "GCA").split(".")[0]
    if not base_api_url.endswith("/"):
        base_api_url += "/"
    response = requests.get(f"{base_api_url}{gb_accession}", timeout=60)
    response.raise_for_status()
    entry = ElementTree.fromstring(response.text)
    taxon_node = entry.find(".//TAXON")
    if taxon_node is None:
        raise MissingNodeError("Cannot find the TAXON node")
    # Fetch taxon ID, scientific_name and strain
    taxon_id = _get_node_text(taxon_node, "TAXON_ID")
    scientific_name = _get_node_text(taxon_node, "SCIENTIFIC_NAME")
    strain = _get_node_text(taxon_node, "STRAIN", optional=True)
    taxonomy = {
        # Ignore arg-type check in the following line since taxon_id cannot be None
        "taxon_id": int(taxon_id),  # type: ignore[arg-type]
        "scientific_name": scientific_name,
    }
    if strain:
        taxonomy["strain"] = strain
    return taxonomy


def _get_node_text(node: Optional[Element], tag: str, optional: bool = False) -> Optional[str]:
    """Returns the value of the field matching the provided tag inside `node`.

    If the tag is not present and `optional` is True, returns `None` instead.

    Args:
        node: Node of an XML tree.
        tag: Tag to fetch within the node.
        optional: Do not raise an exception if the tag does not exist.

    Raises:
        MissingNodeError: If no node is provided or if the tag is missing (if `optional == False`).
    """
    if node is None:
        raise MissingNodeError(f"No node provided to look for '{tag}'")
    tag_text = node.findtext(tag)
    if not optional and (tag_text is None):
        raise MissingNodeError(f"No node found for tag '{tag}'")
    return tag_text


def prepare_genome_metadata(
    input_file: PathLike,
    output_file: PathLike,
    gff3_file: Optional[PathLike] = None,
    base_api_url: str = DEFAULT_API_URL,
    mock_run: bool = False,
) -> None:
    """Updates the genome metadata JSON file with additional information.

    In particular, more information is added about the provider, the assembly and its gene build version,
    and the taxonomy.

    Args:
        input_file: Path to JSON file with genome metadata.
        output_file: Output directory where to generate the final `genome.json` file.
        gff3_file: Path to GFF3 file to use as annotation source for this genome.
        base_api_url: Base API URL to fetch the taxonomy data from.
        mock_run: Do not call external taxonomy service.

    """
    genome_data = get_json(input_file)
    # Amend any missing metadata
    add_provider(genome_data, gff3_file)
    add_assembly_version(genome_data)
    add_genebuild_metadata(genome_data)
    if mock_run:
        genome_data["species"].setdefault("taxonomy_id", 9999999)
    else:
        add_species_metadata(genome_data, base_api_url)
    # Dump updated genome metadata
    print_json(output_file, genome_data)


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path("--input_file", required=True, help="Genome metadata JSON file")
    parser.add_argument_dst_path(
        "--output_file", required=True, help="Output path for the new genome metadata file"
    )
    parser.add_argument_src_path("--gff3_file", help="GFF3 file to use as annotation source")
    parser.add_argument(
        "--base_api_url", default=DEFAULT_API_URL, help="API URL to fetch the taxonomy data from"
    )
    parser.add_argument("--mock_run", action="store_true", help="Do not call external APIs")
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    prepare_genome_metadata(
        input_file=args.input_file,
        output_file=args.output_file,
        gff3_file=args.gff3_file,
        base_api_url=args.base_api_url,
        mock_run=args.mock_run,
    )
