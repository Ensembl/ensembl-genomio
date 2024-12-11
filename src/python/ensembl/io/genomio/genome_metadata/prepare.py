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
    "prepare_genome_metadata",
    "PROVIDER_DATA",
    "MissingNodeError",
    "MetadataError",
]

import datetime
from os import PathLike
from typing import Dict

import ensembl.io.genomio
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


class MissingNodeError(Exception):
    """When a taxon XML node cannot be found."""


class MetadataError(Exception):
    """When a metadata value is not expected."""


def add_provider(genome_metadata: Dict, ncbi_data: Dict) -> None:
    """Updates the genome metadata adding provider information for assembly and gene models.

    Assembly provider metadata will only be added if it is missing, i.e. neither `"provider_name"` or
    `"provider_url"` are present. The gene model metadata will only be added if `gff3_file` is provided.

    Args:
        genome_data: Genome information of assembly, accession and annotation.
        ncbi_data: Report data from NCBI datasets.

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
    if ("provider_name" not in assembly) and ("provider_url" not in assembly):
        assembly["provider_name"] = provider["assembly"]["provider_name"]
        assembly["provider_url"] = f'{provider["assembly"]["provider_url"]}/{accession}'

    # Add annotation provider if there are gene models
    if "annotation_info" in ncbi_data:
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
    if "version" not in assembly:
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
    if "version" not in genebuild:
        genebuild["version"] = current_date
    if "start_date" not in genebuild:
        genebuild["start_date"] = current_date


def add_species_metadata(genome_metadata: Dict, ncbi_data: Dict) -> None:
    """Adds taxonomy ID, scientific name and strain (if present) from the NCBI dataset report.

    Args:
        genome_metadata: Genome information of assembly, accession and annotation.
        ncbi_data: Report data from NCBI datasets.

    """
    species = genome_metadata.setdefault("species", {})
    try:
        organism = ncbi_data["organism"]
    except KeyError:
        return

    if "tax_id" in organism:
        species.setdefault("taxonomy_id", organism["tax_id"])
    if "organism_name" in organism:
        species.setdefault("scientific_name", organism["organism_name"])

    try:
        species.setdefault("strain", organism["infraspecific_names"]["strain"])
    except KeyError:
        pass


def prepare_genome_metadata(
    input_file: PathLike,
    output_file: PathLike,
    ncbi_meta: PathLike,
) -> None:
    """Updates the genome metadata JSON file with additional information.

    In particular, more information is added about the provider, the assembly and its gene build version,
    and the taxonomy.

    Args:
        input_file: Path to JSON file with genome metadata.
        output_file: Output directory where to generate the final `genome.json` file.
        ncbi_meta: JSON file from NCBI datasets.

    """
    genome_data = get_json(input_file)
    ncbi_data = {}
    if ncbi_meta:
        ncbi_data = get_json(ncbi_meta)["reports"][0]

    # Amend any missing metadata
    add_provider(genome_data, ncbi_data)
    add_assembly_version(genome_data)
    add_genebuild_metadata(genome_data)
    add_species_metadata(genome_data, ncbi_data)
    # Dump updated genome metadata
    print_json(output_file, genome_data)


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path("--input_file", required=True, help="Genome metadata JSON file")
    parser.add_argument_dst_path(
        "--output_file", required=True, help="Output path for the new genome metadata file"
    )
    parser.add_argument_src_path(
        "--ncbi_meta", required=True, help="JSON file from NCBI datasets for this genome."
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    prepare_genome_metadata(
        input_file=args.input_file, output_file=args.output_file, ncbi_meta=args.ncbi_meta
    )
