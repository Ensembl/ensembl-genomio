#!/usr/bin/env python3
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

import argparse
import re
from typing import Dict, List, Any

import mysql.connector


class CoreServer:
    """Basic interface to a MySQL server with core databases.

    Allows to get a list of databases and provides access to them.

    Attributes:
        host
        port
        user
        password (optional)

    To connect to a specific database:
    1) Create the core server object
    2) Set the database with core_server.set_database("dbname")
    3) Retrieve a cursor with core_server.get_cursor()
    """

    def __init__(self, host: str, port: str, user: str, password: str = "") -> None:
        self.host = host
        self.port = port
        self.user = user
        self.password = password
        self._connector: Any = None

        # Start a connection directly
        self.connect()

    def connect(self) -> None:
        """Create a connection to the database."""
        self._connector = mysql.connector.connect(
            user=self.user, passwd=self.password, host=self.host, port=self.port
        )

    def set_database(self, db_name: str) -> None:
        self._connector.database = db_name

    def get_cursor(self):
        return self._connector.cursor()

    def get_all_cores(self) -> List[str]:
        """Query the server and retrieve all databases that look like Ensembl cores."""

        query = "SHOW DATABASES LIKE '%_core_%'"

        cursor = self.get_cursor()
        cursor.execute(query)

        dbs = []
        for db in cursor:
            dbs.append(db[0])
        return dbs

    def get_cores(self, prefix: str = "", build: str = "", version: str = "", dbname_re: str = "") -> List[str]:
        """Provide a list of core databases, filtered if requested.
        Args:
            prefix: filter by prefix (no _ is added automatically)
            build: filter by build
            version: filter by Ensembl version
            dbname_re: filter by dbname regular expression

        Returns:
            A list of database names
        """
        dbs = []

        dbs = self.get_all_cores()

        if prefix:
            dbs = [db for db in dbs if db.startswith(f"{prefix}")]
        if dbname_re:
            dbname_m = re.compile(dbname_re)
            dbs = list(filter(dbname_m.search, dbs))
        if build:
            dbs = [db for db in dbs if re.search(rf"_core_{build}_\d+_\d+$", db)]
        if version:
            dbs = [db for db in dbs if re.search(rf"_core_\d+_{version}_\d+$", db)]

        return dbs

    def get_db_metadata(self) -> Dict[str, List]:
        """Retrieve all metadata from a database.

        Returns:
            A dict of with key meta_key, and value=List of meta_value.

        """
        query = "SELECT meta_key, meta_value FROM meta"
        cursor = self.get_cursor()
        cursor.execute(query)

        metadata: Dict[str, List] = {}
        for row in cursor:
            meta_key, meta_value = row
            if meta_key in metadata:
                metadata[meta_key].append(meta_value)
            else:
                metadata[meta_key] = [meta_value]

        return metadata

    def get_table_data(self, table: str, fields: List[str], constraints: str = "") -> List[Dict]:
        """Retrieve all rows from a table.

        Returns:
            A List containing a dict for each row.

        """
        fields_str = ", ".join(fields)
        query = f"SELECT {fields_str} FROM {table}"
        if constraints:
            query = f"{query} WHERE {constraints}"
        cursor = self.get_cursor()
        cursor.execute(query)

        rows = []
        for row in cursor:
            row_data = dict(zip(fields, list(row)))
            rows.append(row_data)

        return rows


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get a list of Ensembl databases")

    parser.add_argument("--host", type=str, required=True, help="Server hostname")
    parser.add_argument("--port", type=str, required=True, help="Server port")
    parser.add_argument("--user", type=str, required=True, help="Server user")
    parser.add_argument("--password", type=str, help="Server password")

    # Optional
    parser.add_argument("--prefix", type=str, help="Prefix for the databases")
    parser.add_argument("--build", type=str, help="Build of the databases")
    parser.add_argument("--version", type=str, help="Ensembl version of the databases")
    args = parser.parse_args()

    # Start
    factory = CoreServer(host=args.host, port=args.port, user=args.user, password=args.password)
    dbs = factory.get_cores(prefix=args.prefix, build=args.build, version=args.version)
    print("\n".join(dbs))
