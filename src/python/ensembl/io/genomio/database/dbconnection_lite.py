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

import logging
import re
from typing import Any, Dict, List, Optional

from sqlalchemy import select
from sqlalchemy.orm import Session

from ensembl.core.models import Meta
from ensembl.utils.database import DBConnection, StrURL

_DB_PATTERN_RELEASE = re.compile(r".+_(?:core|otherfeatures|variation)_(?P<release>\d+)_\d+_\d+")


class DBConnectionLite(DBConnection):
    """Extension to get metadata directly from a database, assuming it has a metadata table."""

    def __init__(self, url: StrURL, reflect: bool = False, **kwargs: Any) -> None:
        super().__init__(url, reflect, **kwargs)
        self._metadata: Dict[str, List] = {}

    def get_metadata(self) -> Dict[str, List]:
        """Retrieves all metadata from the `meta` table in the database.

        Returns:
            A dict of with key meta_key, and value=List of meta_value.

        """
        self._load_metadata()
        return self._metadata

    def _load_metadata(self) -> None:
        """Caches the metadata values."""

        if self._metadata:
            return

        with Session(self._engine) as session:
            meta_stmt = select(Meta)

            for meta_row in session.scalars(meta_stmt).unique().all():
                meta_key = meta_row.meta_key
                meta_value = meta_row.meta_value
                if meta_key in self._metadata:
                    self._metadata[meta_key].append(meta_value)
                else:
                    self._metadata[meta_key] = [meta_value]

    def get_meta_value(self, meta_key: str) -> Optional[str]:
        """Returns the first meta_value for a given meta_key."""

        self._load_metadata()
        try:
            return self._metadata[meta_key][0]
        except KeyError:
            logging.debug(f"No meta_key {meta_key}")
            return None

    def get_project_release(self) -> str:
        """Returns the project release number from the database name. Returns empty string if not found."""

        match = re.search(_DB_PATTERN_RELEASE, self.db_name)
        if match:
            return match.group(1)
        return ""
