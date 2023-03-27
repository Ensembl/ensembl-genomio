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
"""TODO"""

__all__ = ["add_provider", "add_assembly_version", "add_genebuild_metadata", "add_species_metadata",
           "get_taxonomy_from_accession", "prepare_genome_metadata", "PROVIDER_DATA", "DEFAULT_API_URL"]

import datetime
import os
from os import PathLike
from pathlib import Path
import requests
from typing import Dict, Optional
from xml.etree import ElementTree
from xml.etree.ElementTree import Element

import argschema

from ensembl.io.genomio.utils import get_json, print_json


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


def add_provider(genome_data: Dict, gff3_raw: Optional[PathLike] = None) -> None:
    """Adds provider metadata for assembly and gene models in `genome_data`.

    Assembly provider metadata will only be added if it is missing, i.e. neither ``provider_name`` or
    ``provider_url`` are present. The gene model metadata will only be added if `gff3_raw` is provided.

    Args:
        genome_data: Genome information of assembly, accession and annotation.
        gff3_raw: 

    """
    # Get accession provider
    accession = genome_data["assembly"]["accession"]
    if accession.startswith("GCF"):
        provider = PROVIDER_DATA["RefSeq"]
    elif accession.startswith("GCA"):
        provider = PROVIDER_DATA["GenBank"]
    else:
        raise Exception(f"Accession doesn't look like an INSDC or RefSeq accession: {accession}")
    # Add assembly provider (if missing)
    assembly = genome_data["assembly"]
    if (not "provider_name" in assembly) and (not "provider_url" in assembly):
        assembly["provider_name"] = provider["assembly"]["provider_name"]
        assembly["provider_url"] = provider["assembly"]["provider_url"]
    # Add annotation provider if there are gene models
    if gff3_raw:
        annotation = {}
        if "annotation" in genome_data:
            annotation = genome_data["annotation"]
        if (not "provider_name" in annotation) and (not "provider_url" in annotation):
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
    assembly = genome_data["assembly"]
    genebuild = genome_data["genebuild"]
    current_date = datetime.date.today().isoformat()
    if not "version" in genebuild:
        genebuild["version"] = current_date
    if not "start_date" in genebuild:
        genebuild["start_date"] = current_date


def add_species_metadata(genome_data: Dict, base_api_url: Optional[str] = DEFAULT_API_URL) -> None:
    """Adds missing species metadata based on the genome's accession.

    The ``taxonomy_id``, ``strain`` and ``scientific_name`` will be fetched from the taxonomy information
    linked to the given accession.
    
    Args:
        genome_data: Genome information of assembly, accession and annotation.
        base_api_url: Base API URL to fetch the accession's taxonomy data from.

    """
    species = genome_data["species"]
    if not "taxonomy_id" in species:
        accession = genome_data["assembly"]["accession"]
        taxonomy = get_taxonomy_from_accession(accession, base_api_url)
        species["taxonomy_id"] = taxonomy["taxon_id"]
        if (not "strain" in species) and ("strain" in taxonomy):
            species["strain"] = taxonomy["strain"]
        if not "scientific_name" in species:
            species["scientific_name"] = taxonomy["scientific_name"]


def get_taxonomy_from_accession(accession: str, base_api_url: Optional[str] = DEFAULT_API_URL) -> Dict:
    """Returns the taxonomy metadata associated to the given accession.
    
    Args:
        accession: INSDC accession ID.
        base_api_url: Base API URL to fetch the accession's taxonomy data from.

    Returns:
        Dictionary with key-value pairs for ``taxon_id`` and ``scientific_name``. ``strain`` will be added
        only if present in the fetched taxonomy data.

    Raises:
        Exception: If ``TAXON_ID`` or ``SCIENTIFIC_NAME`` are missing in the taxonomy data fetched.

    """
    # Use the GenBank accession without version
    gb_accession = accession.replace("GCF", "GCA").split(".")[0]
    response = requests.get(f"{base_api_url}/{gb_accession}")
    entry = ElementTree.fromstring(response.text)
    taxon_node = entry.find(".//TAXON")
    # Fetch taxon ID, scientific name and strain information
    taxon_id = _get_node_text(taxon_node, "TAXON_ID")
    if not taxon_id:
        raise Exception(f"No 'TAXON_ID' found for accession {accession}")
    scientific_name = _get_node_text(taxon_node, "SCIENTIFIC_NAME")
    if not scientific_name:
        raise Exception(f"No 'SCIENTIFIC_NAME' found for accession {accession}")
    strain = _get_node_text(taxon_node, "STRAIN")

    taxonomy = {
        "taxon_id": int(taxon_id),
        "scientific_name": scientific_name,
    }
    if strain:
        taxonomy["strain"] = strain
    return taxonomy


def _get_node_text(node: Element, tag: str) -> Optional[str]:
    """Returns the value of the field matching the provided tag inside `node`.

    Args:
        node: Node of an XML tree.
        tag: Tag to fetch within the node.

    """
    try:
        return node.find(tag).text
    except:
        return


def prepare_genome_metadata(
    json_path: PathLike,
    output_dir: Optional[PathLike] = "",
    gff3_raw: Optional[PathLike] = None,
    base_api_url: Optional[str] = DEFAULT_API_URL,
) -> None:
    """TODO

    Args:
        json_file: Path to JSON file with genome metadata.
        output_dir: 
        gff3_raw: 
        base_api_url: Base API URL to fetch the accession's taxonomy data from.

    """
    genome_data = get_json(json_path)
    # Amend any missing metadata
    add_provider(genome_data, gff3_raw)
    add_assembly_version(genome_data)
    add_genebuild_metadata(genome_data)
    add_species_metadata(genome_data, base_api_url)
    # Create output directory
    # accession = genome_data["assembly"]["accession"]
    accession = genome_data["BRC4"]["organism_abbrev"]
    work_dir = Path(output_dir, accession)
    work_dir.mkdir(parents=True, exist_ok=True)
    # Dump updated genome metadata
    output_path = work_dir / "genome.json"
    print_json(output_path, genome_data)


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by the entry point of this module."""
    json_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Genome metadata JSON file path"}
    )
    output_dir = argschema.fields.OutputDir(
        required=False, dump_default=".",
        metadata={"description": "Folder where to store the updated genome metadata ($PWD by default)"}
    )


def main() -> None:
    """Module's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    prepare_genome_metadata(mod.args["json_file"], mod.args["output_dir"])
