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
"""Interface to a Mysql server with core databases."""

__all__ = ["CoreServer"]

import re
from typing import List, Optional
import logging

import sqlalchemy
from sqlalchemy.engine import URL
from sqlalchemy import text


class CoreServer:
    """Basic interface to a MySQL server with core databases.

    Allows to get a filtered list of databases.
    """

    def __init__(self, server_url: URL) -> None:
        logging.debug(f"Connect to {server_url}")
        self.engine = sqlalchemy.create_engine(server_url)

    def get_all_core_names(self) -> List[str]:
        """Query the server and retrieve all database names that look like Ensembl cores."""

        with self.engine.connect() as connection:
            all_query = connection.execute(text(r"SHOW DATABASES LIKE '%%_core_%%'"))
            dbs = [row[0] for row in all_query.fetchall()]
        logging.info(f"{len(dbs)} core databases on the server")
        return dbs

    def get_cores(
        self,
        *,
        prefix: str = "",
        build: Optional[int] = None,
        version: Optional[int] = None,
        dbname_re: str = "",
        db_list: Optional[List[str]] = None,
    ) -> List[str]:
        """Returns a list of core databases, filtered if requested.

        Args:
            prefix: Filter by prefix (no "_" is added automatically).
            build: Filter by VEuPathDB build number.
            version: Filter by Ensembl version.
            dbname_re: Filter by dbname regular expression.
            db_list: Explicit list of database names.
        """
        dbs = []

        dbs = self.get_all_core_names()

        # Check if there are databases returned from query to host
        if not dbs:
            logging.warning("No databases returned from query")

        if db_list:
            logging.debug(f"Filter with db list: {db_list}")
            dbs = [db for db in dbs if db in db_list]
        if prefix:
            dbs = [db for db in dbs if db.startswith(f"{prefix}")]
        if dbname_re:
            dbname_m = re.compile(dbname_re)
            dbs = list(filter(dbname_m.search, dbs))
        if build is not None:
            dbs = [db for db in dbs if re.search(rf"_core_{build}_\d+_\d+$", db)]
        if version is not None:
            dbs = [db for db in dbs if re.search(rf"_core_\d+_{version}_\d+$", db)]

        logging.info(f"{len(dbs)} core databases remain after filtering")

        return dbs
