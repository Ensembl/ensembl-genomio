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
"""Expand the genome_metadata with more details for:
the provider, assembly and gene build version, and the taxonomy.
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
import logging
from os import PathLike
from typing import Dict, Optional
from xml.etree import ElementTree
from xml.etree.ElementTree import Element

import requests
from requests.exceptions import ReadTimeout

from ensembl.io.genomio.utils import get_json, print_json
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


PROVIDER_DATA = {
    "GenBank": {
        "assembly": {
            "provider_name": "GenBank",
            "provider_url": "https://www.ncbi.nlm.nih.gov/assembly",
        },
        "annotation": {
            "provider_name": "GenBank",
            "provider_url": "https://www.ncbi.nlm.nih.gov/assembly",
        },
    },
    "RefSeq": {
        "assembly": {
            "provider_name": "RefSeq",
            "provider_url": "https://www.ncbi.nlm.nih.gov/refseq",
        },
        "annotation": {
            "provider_name": "RefSeq",
            "provider_url": "https://www.ncbi.nlm.nih.gov/refseq",
        },
    },
}
DEFAULT_API_URL = "https://www.ebi.ac.uk/ena/browser/api/xml/"


class MissingNodeError(Exception):
    """When a taxon XML node cannot be found."""


class MetadataError(Exception):
    """When a metadata value is not expected."""


def add_provider(genome_data: Dict, gff3_file: Optional[PathLike] = None) -> None:
    """Adds provider metadata for assembly and gene models in `genome_data`.

    Assembly provider metadata will only be added if it is missing, i.e. neither ``provider_name`` or
    ``provider_url`` are present. The gene model metadata will only be added if `gff3_file` is provided.

    Args:
        genome_data: Genome information of assembly, accession and annotation.
        gff3_file: Path to GFF3 file to use as annotation source for this genome.

    """
    # Get accession provider
    accession = genome_data["assembly"]["accession"]
    if accession.startswith("GCF"):
        provider = PROVIDER_DATA["RefSeq"]
    elif accession.startswith("GCA"):
        provider = PROVIDER_DATA["GenBank"]
    else:
        raise MetadataError(f"Accession doesn't look like an INSDC or RefSeq accession: {accession}")

    # Add assembly provider (if missing)
    assembly = genome_data["assembly"]
    if (not "provider_name" in assembly) and (not "provider_url" in assembly):
        assembly["provider_name"] = provider["assembly"]["provider_name"]
        assembly["provider_url"] = provider["assembly"]["provider_url"]

    # Add annotation provider if there are gene models
    if gff3_file:
        annotation = {}
        if "annotation" in genome_data:
            annotation = genome_data["annotation"]
        if ("provider_name" not in annotation) and ("provider_url" not in annotation):
            annotation["provider_name"] = provider["annotation"]["provider_name"]
            annotation["provider_url"] = provider["annotation"]["provider_url"]
        genome_data["annotation"] = annotation


def add_assembly_version(genome_data: Dict) -> None:
    """Adds version number to the genome's assembly if one is not present already.

    Args:
        genome_data: Genome information of assembly, accession and annotation.

    """
    assembly = genome_data["assembly"]
    if not "version" in assembly:
        accession = assembly["accession"]
        values = accession.split(".")
        if (len(values) == 2) and values[1]:
            assembly["version"] = int(values[1])


def add_genebuild_metadata(genome_data: Dict) -> None:
    """Adds missing genebuild metadata.

    The default convention is to use the current date as ``version`` and ``start_date``.

    Args:
        genome_data: Genome information of assembly, accession and annotation.

    """
    genebuild = genome_data["genebuild"]
    current_date = datetime.date.today().isoformat()
    if not "version" in genebuild:
        genebuild["version"] = current_date
    if not "start_date" in genebuild:
        genebuild["start_date"] = current_date


def add_species_metadata(genome_data: Dict, base_api_url: str = DEFAULT_API_URL) -> None:
    """Adds missing species metadata based on the genome's accession.

    The ``taxonomy_id``, ``strain`` and ``scientific_name`` will be fetched from the taxonomy information
    linked to the given accession.

    Args:
        genome_data: Genome information of assembly, accession and annotation.
        base_api_url: Base API URL to fetch the taxonomy data from.

    """
    species = genome_data["species"]
    if not "taxonomy_id" in species:
        try:
            accession = genome_data["assembly"]["accession"]
            taxonomy = get_taxonomy_from_accession(accession, base_api_url)
            species["taxonomy_id"] = taxonomy["taxon_id"]
            if (not "strain" in species) and ("strain" in taxonomy):
                species["strain"] = taxonomy["strain"]
            if not "scientific_name" in species:
                species["scientific_name"] = taxonomy["scientific_name"]
        except KeyError:
            logging.warning("Could not extract taxonomy data")



def get_taxonomy_from_accession(accession: str, base_api_url: str = DEFAULT_API_URL) -> Dict:
    """Returns the taxonomy metadata associated to the given accession.

    Args:
        accession: INSDC accession ID.
        base_api_url: Base API URL to fetch the taxonomy data from.

    Returns:
        Dictionary with key-value pairs for ``taxon_id`` and ``scientific_name``. ``strain`` will be added
        only if present in the fetched taxonomy data.

    Raises:
        MissinDataException: If ``TAXON_ID`` or ``SCIENTIFIC_NAME`` are missing in the taxonomy data fetched.

    """
    # Use the GenBank accession without version
    gb_accession = accession.replace("GCF", "GCA").split(".")[0]

    try:
        response = requests.get(f"{base_api_url}/{gb_accession}", timeout=5)
        response.raise_for_status()

        entry = ElementTree.fromstring(response.text)

        taxon_node = entry.find(".//TAXON")
        if taxon_node is None:
            raise MissingNodeError("Can't find the TAXON node")

        # Fetch taxon ID, scientific_name and strain
        taxon_id = _get_node_text(taxon_node, "TAXON_ID")
        scientific_name = _get_node_text(taxon_node, "SCIENTIFIC_NAME")
        strain = _get_node_text(taxon_node, "STRAIN", optional=True)

        if taxon_id and scientific_name:
            taxonomy = {
                "taxon_id": int(taxon_id),
                "scientific_name": scientific_name,
            }
        if strain:
            taxonomy["strain"] = strain
    except ReadTimeout:
        logging.warning("Can't get taxonomy from ENA")
        taxonomy = {"taxon_id": 0, "scientific_name": ""}

    return taxonomy


def _get_node_text(node: Element, tag: str, optional: bool = False) -> Optional[str]:
    """Returns the value of the field matching the provided tag inside `node`.
    By default raise a MissingNodeException if the tag is not found.
    If optional is True and no tag is found, return None.

    Args:
        node: Node of an XML tree.
        tag: Tag to fetch within the node.
        optional: Don't raise an exception if the tag doesn't exist.

    """
    if node is None:
        raise MissingNodeError(f"No node provided to look for {tag}")
    tag_node = node.find(tag)

    if tag_node is not None:
        return tag_node.text
    if optional:
        return None
    raise MissingNodeError(f"No node found for tag {tag}")


def prepare_genome_metadata(
    input_file: PathLike,
    output_file: PathLike,
    gff3_file: Optional[PathLike] = None,
    base_api_url: str = DEFAULT_API_URL,
) -> None:
    """Updates the genome metadata JSON file with additional information.

    In particular, more information is added about the provider, the assembly and its gene build version,
    and the taxonomy.

    Args:
        input_file: Path to JSON file with genome metadata.
        output_file: Output directory where to generate the final `genome.json` file.
        gff3_file: Path to GFF3 file to use as annotation source for this genome.
        base_api_url: Base API URL to fetch the taxonomy data from.

    """
    genome_data = get_json(input_file)
    # Amend any missing metadata
    add_provider(genome_data, gff3_file)
    add_assembly_version(genome_data)
    add_genebuild_metadata(genome_data)
    add_species_metadata(genome_data, base_api_url)
    # Dump updated genome metadata
    print_json(output_file, genome_data)


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(
        description=(
            "Add information about provider, taxonomy and assembly and gene build version to the genome "
            "metadata file."
        )
    )
    parser.add_argument_src_path("--input_file", required=True, help="Genome metadata JSON file")
    parser.add_argument_dst_path(
        "--output_file", required=True, help="Output path for the new genome metadata file"
    )
    parser.add_argument_src_path("--gff3_file", help="GFF3 file to use as annotation source")
    parser.add_argument(
        "--base_api_url", default=DEFAULT_API_URL, help="API URL to fetch the taxonomy data from"
    )
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    prepare_genome_metadata(
        input_file=args.input_file,
        output_file=args.output_file,
        gff3_file=args.gff3_file,
        base_api_url=args.base_api_url,
    )
