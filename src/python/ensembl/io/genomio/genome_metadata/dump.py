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
"""Generates a JSON file representing the genome metadata from a core database."""

__all__ = [
    "get_genome_metadata",
    "filter_genome_meta",
    "check_assembly_version",
    "check_genebuild_version",
]

import json
from typing import Any, Dict, Type
import logging

from sqlalchemy import select
from sqlalchemy.orm import Session

from ensembl.core.models import Meta
from ensembl.io.genomio.database import DBConnectionLite
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


METADATA_FILTER: Dict[str, Dict[str, Type]] = {
    "added_seq": {"region_name": str},
    "annotation": {"provider_name": str, "provider_url": str},
    "assembly": {
        "accession": str,
        "date": str,
        "name": str,
        "provider_name": str,
        "provider_url": str,
        "version": int,
    },
    "BRC4": {"organism_abbrev": str, "component": str},
    "genebuild": {"id": str, "method": str, "method_display": str, "start_date": str, "version": str},
    "species": {
        "alias": str,
        "annotation_source": str,
        "display_name": str,
        "division": str,
        "production_name": str,
        "scientific_name": str,
        "strain": str,
        "taxonomy_id": int,
    },
}


def get_genome_metadata(session: Session) -> Dict[str, Any]:
    """Returns the meta table content from the core database in a nested dictionary.

    Args:
        session: Session for the current core.

    """
    genome_metadata: Dict[str, Any] = {}
    meta_statement = select(Meta)
    for row in session.execute(meta_statement).unique().all():
        meta_key = row[0].meta_key
        meta_value = row[0].meta_value
        (main_key, _, subkey) = meta_key.partition(".")
        # Use empty string as subkey when no "." found to simplify dictionary creation
        if main_key in genome_metadata:
            if subkey in genome_metadata[main_key]:
                genome_metadata[main_key][subkey].append(meta_value)
            else:
                genome_metadata[main_key][subkey] = [meta_value]
        else:
            genome_metadata[main_key] = {subkey: [meta_value]}
    # Parse genome metadata to simplify dictionary and check data consistency
    for main_key, subkeys_dict in genome_metadata.items():
        # Replace single-value lists by the value itself
        for subkey, value in subkeys_dict.items():
            if len(value) == 1:
                subkeys_dict[subkey] = value[0]
        # Remove nested dictionary if it only has "" as key, passing its value to the main key
        if "" in subkeys_dict:
            if len(subkeys_dict) == 1:
                genome_metadata[main_key] = subkeys_dict.pop("")
            else:
                raise ValueError(f"Unexpected meta keys for '{main_key}': {', '.join(subkeys_dict.keys())}")
    return genome_metadata


def filter_genome_meta(genome_metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Returns a filtered metadata dictionary with only the predefined keys in METADATA_FILTER.

    Also converts to expected data types (to follow the genome JSON schema).

    Args:
        genome_metadata: Nested metadata key values from the core metadata table.
    """
    filtered_metadata: Dict[str, Any] = {}
    for key, subfilter in METADATA_FILTER.items():
        if key in genome_metadata:
            filtered_metadata[key] = {}
            for subkey, value_type in subfilter.items():
                if subkey in genome_metadata[key]:
                    value = genome_metadata[key][subkey]
                    if isinstance(value, list):
                        value = [value_type(x) for x in value]
                    else:
                        value = value_type(value)
                    filtered_metadata[key][subkey] = value
    # Check assembly and genebuild versions
    check_assembly_refseq(filtered_metadata)
    check_assembly_version(filtered_metadata)
    check_genebuild_version(filtered_metadata)
    return filtered_metadata


def check_assembly_refseq(gmeta_out: Dict[str, Any]) -> None:
    """Update the GCA accession to use GCF if it is from RefSeq.

    Args:
        genome_metadata: Nested metadata key values from the core metadata table.
    """
    assembly = gmeta_out.get("assembly", {})
    if assembly.get("provider_name", "") == "RefSeq":
        assembly["accession"] = assembly["accession"].replace("GCA", "GCF")


def check_assembly_version(genome_metadata: Dict[str, Any]) -> None:
    """Updates the assembly version of the genome metadata provided.

    If `version` meta key is not and integer or it is not available, the assembly accession's version
    will be used instead.

    Args:
        genome_metadata: Nested metadata key values from the core metadata table.

    Raises:
        ValueError: If both `version` and the assembly accession's version are not integers or are missing.
    """
    assembly = genome_metadata["assembly"]
    version = assembly.get("version")
    # Check the version is an integer
    try:
        assembly["version"] = int(version)
    except (ValueError, TypeError) as exc:
        # Get the version from the assembly accession
        accession = assembly["accession"]
        version = accession.partition(".")[2]
        try:
            assembly["version"] = int(version)
        except ValueError:
            raise ValueError(f"Assembly version is not an integer in {assembly}") from exc
        logging.info(f"Assembly version [v{version}] obtained from assembly accession ({accession}).")
    else:
        logging.info(f'Located version [v{assembly["version"]}] info from meta data.')


def check_genebuild_version(genome_metadata: Dict[str, Any]) -> None:
    """Updates the genebuild version (if not present) from the genebuild ID, removing the latter.

    Args:
        genome_metadata: Nested metadata key values from the core metadata table.

    Raises:
        ValueError: If there is no genebuild version or ID available.
    """
    try:
        genebuild = genome_metadata["genebuild"]
    except KeyError:
        return
    if "version" not in genebuild:
        try:
            genebuild_id = genebuild["id"]
        except KeyError:
            # pylint: disable=raise-missing-from
            raise ValueError("No genebuild version or ID found")
        genome_metadata["genebuild"]["version"] = str(genebuild_id)
    # Drop genebuild ID since there is a genebuild version
    genome_metadata["genebuild"].pop("id", None)


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Fetch the genome metadata from a core database and print it in JSON format."
    )
    parser.add_server_arguments(include_database=True)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnectionLite(args.url)
    with dbc.session_scope() as session:
        genome_meta = get_genome_metadata(session)
        genome_meta = filter_genome_meta(genome_meta)

    print(json.dumps(genome_meta, indent=2, sort_keys=True))
