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
"""Generates one JSON file per metadata type inside `manifest`, including the manifest itself."""

__all__ = ["format_db_data", "get_core_dbs_metadata"]

import argparse
import json
import logging
from pathlib import Path

from sqlalchemy.engine import URL

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args
from .core_server import CoreServer
from .dbconnection_lite import DBConnectionLite


def format_db_data(server_url: URL, dbs: list[str], brc_mode: bool = False) -> list[dict]:
    """Returns a metadata list from the given databases on a server.

    Args:
        server_url: Server URL where all the databases are hosted.
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
        core_db = DBConnectionLite(db_url)

        prod_name = core_db.get_meta_value("species.production_name")
        species = prod_name
        division = core_db.get_meta_value("species.division")
        accession = core_db.get_meta_value("assembly.accession")
        project_release = core_db.get_project_release()

        if brc_mode:
            brc_organism = core_db.get_meta_value("BRC4.organism_abbrev")
            brc_component = core_db.get_meta_value("BRC4.component")
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


def get_core_dbs_metadata(
    server_url: URL,
    *,
    prefix: str = "",
    build: int | None = None,
    version: int | None = None,
    db_regex: str = "",
    db_list: Path | None = None,
    brc_mode: bool = False,
) -> list[dict]:
    """Returns all the metadata fetched for the selected core databases.

    Args:
        server_url: Server URL where the core databases are stored.
        prefix: Filter by prefix (no "_" is added automatically).
        build: Filter by VEuPathDB build number.
        version: Filter by Ensembl version.
        db_regex: Filter by dbname regular expression.
        db_list: Explicit list of database names.
        brc_mode: Enable BRC mode.

    Returns:
        List of dictionaries with 3 keys: "database", "species" and "division".
    """
    db_list_file = None
    if db_list:
        with db_list.open("r") as infile_fh:
            db_list_file = [line.strip() for line in infile_fh]
    # Get all database names
    server = CoreServer(server_url)
    logging.debug("Fetching databases...")
    databases = server.get_cores(
        prefix=prefix, build=build, version=version, dbname_re=db_regex, db_list=db_list_file
    )
    logging.info(f"Got {len(databases)} databases")
    logging.debug("\n".join(databases))
    return format_db_data(server_url, databases, brc_mode)


def parse_args(arg_list: list[str] | None) -> argparse.Namespace:
    """Return a populated namespace with the arguments parsed from a list or from the command line.

    Args:
        arg_list: List of arguments to parse. If `None`, grab them from the command line.

    """
    parser = ArgumentParser(description=__doc__)
    parser.add_server_arguments()
    # Add filter arguments
    parser.add_argument("--prefix", default="", help="Prefix to filter the databases")
    parser.add_argument("--build", type=int, default=None, help="Build to filter the databases")
    parser.add_argument("--release", type=int, default=None, help="EnsEMBL release to filter the databases")
    parser.add_argument("--db_regex", default="", help="Regular expression to match database names against")
    parser.add_argument_src_path("--db_list", help="File with one database per line to load")
    # Add flags
    parser.add_argument(
        "--brc_mode",
        action="store_true",
        help="Enable BRC mode, i.e. use organism_abbrev for species, component for division",
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    return parser.parse_args(arg_list)


def main(arg_list: list[str] | None = None) -> None:
    """Main script entry-point.

    Args:
        arg_list: Arguments to parse passing list to parse_args().

    """
    args = parse_args(arg_list)
    init_logging_with_args(args)

    databases_data = get_core_dbs_metadata(
        server_url=args.url,
        prefix=args.prefix,
        build=args.build,
        version=args.release,
        db_regex=args.db_regex,
        db_list=args.db_list,
        brc_mode=args.brc_mode,
    )
    print(json.dumps(databases_data, sort_keys=True, indent=4))
