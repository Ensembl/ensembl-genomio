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

__all__ = ["format_db_data", "get_metadata_value"]

import json
from os import PathLike
from pathlib import Path
import re
from typing import Dict, List, Optional
import logging

from sqlalchemy.engine import URL

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args
from .core_server import CoreServer
from .core_database import CoreDatabase


_DB_PATTERN = re.compile(r".+_core_(\d+)_\d+_\d+")


def format_db_data(server_url: URL, dbs: List[str], brc_mode: bool = False) -> List[Dict]:
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
    for db_name in dbs:
        logging.debug(f"Get metadata for {db_name}")
        db_url = server_url.set(database=db_name)
        core_db = CoreDatabase(db_url)
        metadata = core_db.get_metadata()

        prod_name = get_metadata_value(metadata, "species.production_name")
        species = prod_name
        division = get_metadata_value(metadata, "species.division")
        accession = get_metadata_value(metadata, "assembly.accession")
        project_release = _get_project_release(db_name)

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
            "host": db_url.host,
            "user": db_url.username,
            "port": db_url.port,
            "password": db_url.password,
            "database": db_url.database,
        }
        db_data = {
            "server": server_data,
            "production_name": prod_name,
            "species": species,
            "division": division,
            "accession": accession,
            "release": project_release,
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

    match = re.search(_DB_PATTERN, db_name)
    if match:
        return match.group(1)
    return ""


def _load_multine_file(infile: PathLike) -> List[str]:
    data_list = []
    with Path(infile).open("r") as infile_fh:
        data_list = [line.strip() for line in infile_fh]
    return data_list


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
    parser.add_argument_src_path("--db_list", help="File with one database per line to load")
    # Add flags
    parser.add_argument(
        "--brc_mode",
        action="store_true",
        help="Enable BRC mode, i.e. use organism_abbrev for species, component for division",
    )
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    db_list_file = None
    if args.db_list:
        db_list_file = _load_multine_file(args.db_list)

    # Get all db names
    server_url = URL(
        drivername="mysql",
        host=args.host,
        port=args.port,
        username=args.user,
        password=args.password,
    )
    server = CoreServer(server_url)
    logging.debug("Get databases...")
    databases = server.get_cores(
        prefix=args.prefix,
        build=args.build,
        version=args.version,
        dbname_re=args.db_regex,
        db_list=db_list_file,
    )
    logging.info(f"Got {len(databases)} databases")
    logging.debug("\n".join(databases))

    # Get all metadata for those databases
    databases_data = format_db_data(server_url, databases, args.brc_mode)
    print(json.dumps(databases_data, sort_keys=True, indent=4))


if __name__ == "__main__":
    main()
