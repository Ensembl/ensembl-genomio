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
"""CoreDatabase interface to an Ensembl core database."""

__all__ = ["CoreDatabase"]

from typing import Dict, List

from sqlalchemy import create_engine, select
from sqlalchemy.engine import URL

from ensembl.database import DBConnection
from ensembl.core.models import Meta


class CoreDatabase(DBConnection):
    """Add some useful interface for an Ensembl core database."""

    def __init__(self, url: URL, **kwargs) -> None:
        # super().__init__(url, **kwargs)
        self._engine = create_engine(url, **kwargs)
        if self.schema_type != "core":
            raise TypeError(f"This is not a core database ({self.schema_type})")

    @property
    def schema_type(self) -> str:
        with self.session_scope() as session:
            meta_stmt = select(Meta.meta_value).where(Meta.meta_key == "schema_type")
            meta_row = session.execute(meta_stmt).first()
            return meta_row.meta_value

    def get_metadata(self) -> Dict[str, List]:
        """Retrieve all metadata from the `meta` table in the core database.

        Returns:
            A dict of with key meta_key, and value=List of meta_value.

        """
        with self.session_scope() as session:
            meta_stmt = select(Meta)

            metadata: Dict[str, List] = {}
            for meta_row in session.scalars(meta_stmt).unique().all():
                meta_key = meta_row.meta_key
                meta_value = meta_row.meta_value
                if meta_key in metadata:
                    metadata[meta_key].append(meta_value)
                else:
                    metadata[meta_key] = [meta_value]

        return metadata
