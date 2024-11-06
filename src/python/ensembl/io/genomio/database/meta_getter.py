#!/usr/bin/env python
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
"""A simple helper script to connect to a core database and retrieve a single meta_value 
or multiple meta_value and dump meta_key/value pairs to stdout / JSON."""

__all__ = ["get_meta_values"]

import argparse
import logging
import json
from pathlib import PosixPath
from pathlib import Path

from sqlalchemy.engine import URL

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils import StrPath
from ensembl.utils.logging import init_logging_with_args
from .dbconnection_lite import DBConnectionLite


def get_meta_values(server_url: URL, db_name: str, meta_keys: StrPath | list[str]) -> dict[str, str]:
    """Returns a set of meta values based on set of 1 or more input DB meta_keys.

    Args:
        server_url: Server URL where the core databases are stored.
        db_name: Name of the target DB to query.
        query_meta_keys: The meta table 'meta_key' list to query.

    """
    db_url = server_url.set(database=db_name)
    core_db = DBConnectionLite(db_url)
    query_meta_keys = []
    unpopulated_meta_keys = []
    meta_values_located = {}
    input_keys_count = 0
    meta_populated = False

    # Check input type and populated query list
    if isinstance(input_queries, PosixPath):
        with Path(input_queries).open(mode="r", encoding="UTF-8") as fh:
            for line in fh.readlines():
                meta_key = line.strip()
                query_meta_keys.append(meta_key)
    elif isinstance(input_queries, list):
        query_meta_keys = input_queries

    # Loop over input meta_key(s) and query DB
    if len(query_meta_keys) >= 1:
        for meta_key in query_meta_keys:
            input_keys_count += 1
            meta_value = core_db.get_meta_value(f"{meta_key}")

            if meta_value is not None:
                meta_values_located[f"{meta_key}"] = meta_value
            else:
                unpopulated_meta_keys.append(f"{meta_key}")
                logging.info(f"Meta query returned no entry on meta_key: '{meta_key}'")
    else:
        logging.warning(f"No meta_keys found in input file {query_meta_keys}")
        return None

    # Now assess what meta info was recovered and dump to JSON
    total_queries_located = len(meta_values_located.items())
    if total_queries_located == input_keys_count:
        meta_populated = True
    elif (total_queries_located >= 1) and (total_queries_located < input_keys_count):
        meta_populated = True
        logging.info(
            f"Some query meta_keys missing [Queries (input: {input_keys_count}) vs (Located: {total_queries_located})"
        )
        logging.info(f"Missing meta_key(s)-> {unpopulated_meta_keys}")
    else:
        logging.warning("Zero input query meta_keys present/populated.")

    if meta_populated is True:
        meta_values_located["database_name"] = f"{db_name}"
        print(json.dumps(meta_values_located, sort_keys=True, indent="  "))
        return meta_values_located
    else:
        return {}


def parse_args(arg_list: list[str] | None) -> argparse.Namespace:
    """TODO

    Args:
        arg_list: TODO
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_server_arguments()
    parser.add_argument("--database_name", default=None, help="Target database name.")
    parser.add_argument_src_path(
        "--meta_keys_list", help="Input File | List with >=2 meta_keys to query target database."
    )
    parser.add_log_arguments(add_log_file=False)
    return parser.parse_args(arg_list)


def main(arg_list: list[str] | None = None) -> None:
    args = parse_args(arg_list)
    init_logging_with_args(args)

    _ = get_meta_values(server_url=args.url, db_name=args.database_name, input_queries=args.meta_keys_list)
