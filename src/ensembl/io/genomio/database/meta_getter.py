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
"""Connect to a core database and retrieve a meta_key:meta_value pair(s) and dump meta_key/value
pairs to stdout / JSON.
"""

__all__ = ["get_meta_values"]

import argparse
import logging
import json
from pathlib import PosixPath
from pathlib import Path

from sqlalchemy.engine import URL

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils import StrPath
from ensembl.utils.logging import init_logging_with_args
from .dbconnection_lite import DBConnectionLite


def get_meta_values(db_url: URL, meta_keys: StrPath | list[str]) -> dict[str, str]:
    """Returns a set of meta values based on set of 1 or more input DB meta_keys.

    Args:
        db_url: Target core database URL.
        meta_keys: File path with one meta key per line or list of meta keys.

    """
    db_name = db_url.database
    core_db = DBConnectionLite(db_url)
    query_meta_keys = []
    unpopulated_meta_keys = []
    meta_values_located = {}
    input_keys_count = 0
    meta_populated = False

    # Check input type and populated query list
    if isinstance(meta_keys, PosixPath):
        with Path(meta_keys).open(mode="r", encoding="UTF-8") as fh:
            for line in fh.readlines():
                meta_key = line.strip()
                query_meta_keys.append(meta_key)
    elif isinstance(meta_keys, list):
        query_meta_keys = meta_keys

    # Loop over input meta_key(s) and query DB
    for meta_key in query_meta_keys:
        input_keys_count += 1
        meta_value = core_db.get_meta_value(f"{meta_key}")

        if meta_value is not None:
            meta_values_located[f"{meta_key}"] = meta_value
        else:
            unpopulated_meta_keys.append(f"{meta_key}")
            logging.info(f"Meta query returned no entry on meta_key: '{meta_key}'")

    # Now assess what meta info was recovered and dump to JSON
    total_queries_located = len(meta_values_located)
    if total_queries_located >= 1:
        meta_populated = True
        if total_queries_located < input_keys_count:
            logging.info(f"Missing meta_key(s)-> {unpopulated_meta_keys}")
    else:
        logging.warning("Zero input query meta_keys present/populated.")

    if meta_populated:
        meta_values_located["database_name"] = f"{db_name}"
        print(json.dumps(meta_values_located, sort_keys=True, indent=2))
        return meta_values_located
    return {}


def parse_args(arg_list: list[str] | None) -> argparse.Namespace:
    """Return a populated namespace with the arguments parsed from a list or from the command line.

    Args:
        arg_list: List of arguments to parse. If `None`, grab them from the command line.

    """
    parser = ArgumentParser(description=__doc__)
    parser.add_server_arguments(include_database=True, help="server url and core database")
    parser.add_argument_src_path(
        "--meta_keys_list", help="Input File | List with >=2 meta_keys to query target database."
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=False)
    return parser.parse_args(arg_list)


def main(arg_list: list[str] | None = None) -> None:
    """Main script entry-point.

    Args:
        arg_list: Arguments to parse passing list to parse_args().
    """
    args = parse_args(arg_list)
    init_logging_with_args(args)

    _ = get_meta_values(db_url=args.url, meta_keys=args.meta_keys_list)
