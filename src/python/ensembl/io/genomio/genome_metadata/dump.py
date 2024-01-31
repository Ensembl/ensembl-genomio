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

__all__ = ["get_genome_metadata", "filter_genome_meta", "check_assembly_version"]

import json
from typing import Any, Dict
import logging

from sqlalchemy import select
from sqlalchemy.orm import Session

from ensembl.core.models import Meta
from ensembl.database import DBConnection
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def get_genome_metadata(session: Session) -> Dict[str, Any]:
    """Retrieve a select list of metadata from the core database.

    Args:
        session: Session for the current core.

    Returns:
        A nested dict.
    """
    gmeta: Dict[str, Any] = {}

    gmeta_st = select(Meta)
    for row in session.execute(gmeta_st).unique().all():
        dat = row[0]
        meta_key = dat.meta_key
        meta_value = dat.meta_value

        if "." in meta_key:
            (high_key, low_key) = meta_key.split(".")
            if high_key in gmeta:
                if low_key in gmeta[high_key]:
                    gmeta[high_key][low_key].append(meta_value)
                else:
                    gmeta[high_key][low_key] = [meta_value]
            else:
                gmeta[high_key] = {}
                gmeta[high_key][low_key] = [meta_value]
        else:
            if meta_key in gmeta:
                gmeta[meta_key].append(meta_value)
            else:
                gmeta[meta_key] = [meta_value]

    return gmeta


def filter_genome_meta(gmeta: Dict[str, Any]) -> Dict[str, Any]:
    """Returns a filtered metadata dict with only predefined keys.
    Also converts expected numbers to integers (to follow the genome json schema).

    Args:
        gmeta (Dict[str, Any]): Nested metadata key values from the core metadata table.
    """
    meta_list = {
        "species": {
            "taxonomy_id",
            "production_name",
            "scientific_name",
            "strain",
            "display_name",
            "division",
            "alias",
            "annotation_source",
        },
        "assembly": {"accession", "date", "name", "version", "provider_name", "provider_url"},
        "genebuild": {"version", "method", "start_date", "method_display", "id"},
        "annotation": {"provider_name", "provider_url"},
        "BRC4": {"organism_abbrev", "component"},
        "added_seq": {"region_name"},
    }
    is_integer = {"species": {"taxonomy_id"}, "assembly": {"version"}}

    gmeta_out: Dict[str, Any] = {}
    for key1, subkeys in meta_list.items():
        if key1 not in gmeta:
            continue
        if subkeys:
            gmeta_out[key1] = {}
            for key2 in subkeys:
                if key2 not in gmeta[key1]:
                    continue
                value = gmeta[key1][key2]
                if len(value) == 1:
                    value = value[0]
                    if key2 in is_integer.get(key1, {}):
                        value = int(value)
                gmeta_out[key1][key2] = value
        else:
            value = gmeta[key1]
            if len(value) == 1:
                value = value[0]
                if is_integer.get(key1):
                    value = int(value)
            gmeta_out[key1] = value

    check_assembly_version(gmeta_out)
    check_assembly_refseq(gmeta_out)
    check_genebuild_version(gmeta_out)

    return gmeta_out


def check_assembly_refseq(gmeta_out: Dict[str, Any]) -> None:
    """Update the GCA accession to use GCF if it is from RefSeq.

    Args:
        gmeta (Dict[str, Any]): Nested metadata key values from the core metadata table.
    """
    assembly = gmeta_out["assembly"]
    if assembly.get("provider_name", "") == "RefSeq":
        assembly["accession"] = assembly["accession"].replace("GCA", "GCF")


def check_assembly_version(gmeta_out: Dict[str, Any]) -> None:
    """Update the assembly version of the genome metadata provided to use an integer.
    Get the version from the assembly accession as alternative.

    Args:
        gmeta (Dict[str, Any]): Nested metadata key values from the core metadata table.
    """
    assembly = gmeta_out["assembly"]
    version = assembly.get("version")

    # Check the version is an integer
    try:
        assembly["version"] = int(version)
    except (ValueError, TypeError) as exc:
        # Get the version from the assembly accession
        accession = assembly["accession"]
        parts = accession.split(".")
        if len(parts) == 2 and parts[1].isdigit():
            version = parts[1]
            assembly["version"] = int(version)
            logging.info(
                f'Asm version [v{version}] obtained from: assembly accession ({assembly["accession"]}).'
            )
        else:
            raise ValueError(f"Assembly version is not an integer in {assembly}") from exc
    else:
        logging.info(f"Located version [v{int(version)}] info from meta data.")


def check_genebuild_version(metadata: Dict[str, Any]) -> None:
    """Updates the genebuild version (if not present) from the genebuild ID, removing the latter.

    Args:
        metadata: Nested metadata key values from the core metadata table.

    Raises:
        ValueError: If there is no genebuild version or ID available.
    """
    genebuild = metadata.get("genebuild")
    if genebuild is None:
        return
    version = genebuild.get("version")

    # Check there is a version
    if version is None:
        gb_id = genebuild.get("id")
        if gb_id is None:
            raise ValueError("No genebuild version or id")
        metadata["genebuild"]["version"] = str(gb_id)

    if "id" in genebuild:
        del metadata["genebuild"]["id"]


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Fetch the genome metadata from a core database and print it in JSON format."
    )
    parser.add_server_arguments(include_database=True)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnection(args.url)
    with dbc.session_scope() as session:
        genome_meta = get_genome_metadata(session)
        genome_meta = filter_genome_meta(genome_meta)

    print(json.dumps(genome_meta, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
