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
"""Simplified Database interface to an Ensembl database, to make access to metadata easier and faster."""

__all__ = ["DBConnectionLite"]

from typing import Dict, List

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from ensembl.database import DBConnection
from ensembl.core.models import Meta


class DatabaseExtension:
    """Extension to get metadata directly from a database, assuming it has a metadata table. Do use directly.
    """
    def __init__(self, url, **kwargs) -> None:
        self._engine = create_engine(url, **kwargs)
        self._metadata: Dict[str, List] = {}

    def get_metadata(self) -> Dict[str, List]:
        """Retrieves all metadata from the `meta` table in the database.

        Returns:
            A dict of with key meta_key, and value=List of meta_value.

        """
        if self._metadata:
            return self._metadata
        with Session(self._engine) as session:
            meta_stmt = select(Meta)

            for meta_row in session.scalars(meta_stmt).unique().all():
                meta_key = meta_row.meta_key
                meta_value = meta_row.meta_value
                if meta_key in self._metadata:
                    self._metadata[meta_key].append(meta_value)
                else:
                    self._metadata[meta_key] = [meta_value]

        return self._metadata


class DBConnectionLite(DatabaseExtension, DBConnection):
    """DBConnection without the schema loading from DB: faster but some methods will not work."""
