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
"""Generates one JSON file per metadata type inside `manifest`, including the manifest itself.

Can be imported as a module and called as a script as well, with the same parameters and expected outcome.
"""

import json
import re
from typing import Dict, List, Optional

from ensembl.brc4.runnable.core_server import CoreServer
from ensembl.utils.argparse import ArgumentParser


db_pattern = re.compile(r".+_core_(\d+)_(\d+)_\d+")


def format_db_data(server: CoreServer, dbs: List[str], brc_mode: bool = False) -> List[Dict]:
    """Returns metadata from a list of databases on a server.

    Args:
        server: Server where all the databases are hosted.
        dbs: List of database names.
        brc_mode: If true, assign ``BRC4.organism_abbrev`` as the species, and ``BRC4.component`` as the
            division. Otherwise, the species will be ``species.production_name`` and the division will be
            ``species.division``.

    Returns:
        List of dictionaries with 3 keys: "database", "species" and "division".

    """
    databases_data = []
    for db in dbs:
        server.set_database(db)
        metadata = server.get_db_metadata()

        prod_name = get_metadata_value(metadata, "species.production_name")
        species = prod_name
        division = get_metadata_value(metadata, "species.division")
        accession = get_metadata_value(metadata, "assembly.accession")
        project_release = _get_project_release(db)
        ensembl_version = get_metadata_value(metadata, "schema_version")

        if brc_mode:
            brc_organism = get_metadata_value(metadata, "BRC4.organism_abbrev")
            brc_component = get_metadata_value(metadata, "BRC4.component")
            if brc_organism is not None:
                species = brc_organism
            if brc_component is not None:
                division = brc_component

        if not division:
            division = "all"

        server_data = {
            "host": server.host,
            "user": server.user,
            "port": server.port,
            "password": server.password,
            "database": db,
        }
        db_data = {
            "server": server_data,
            "production_name": prod_name,
            "species": species,
            "division": division,
            "accession": accession,
            "release": project_release,
            "ensembl_version": ensembl_version,
        }

        databases_data.append(db_data)
    return databases_data


def get_metadata_value(metadata: Dict[str, List], key: str) -> Optional[str]:
    """Returns the first element in the list assigned to `key` in `metadata`.

    Args:
        metadata: Map of metadata information to lists of values.
        key: Metadata key to search for.

    """
    if (key in metadata) and metadata[key]:
        return metadata[key][0]
    return None


def _get_project_release(db_name: str) -> str:
    """Return the project release number from the database name."""

    match = re.search(db_pattern, db_name)
    if match:
        return match.group(1)
    return ""


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Get the metadata from a list of databases on a server (in JSON format)."
    )
    parser.add_server_arguments()
    # Add filter arguments
    parser.add_argument("--prefix", default="", help="Prefix to filter the databases")
    parser.add_argument("--build", default="", help="Build to filter the databases")
    parser.add_argument("--version", default="", help="EnsEMBL version to filter the databases")
    parser.add_argument("--db_regex", default="", help="Regular expression to match database names against")
    # Add flags
    parser.add_argument(
        "--brc_mode",
        action="store_true",
        help="Enable BRC mode, i.e. use organism_abbrev for species, component for division",
    )
    args = parser.parse_args()

    server = CoreServer(host=args.host, port=args.port, user=args.user, password=args.password)
    databases = server.get_cores(
        prefix=args.prefix, build=args.build, version=args.version, dbname_re=args.db_regex
    )
    databases_data = format_db_data(server, databases, args.brc_mode)
    print(json.dumps(databases_data, sort_keys=True, indent=4))


if __name__ == "__main__":
    main()
