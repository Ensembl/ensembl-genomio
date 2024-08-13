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
"""Get a mapping for external db names."""

__all__ = [
    "DEFAULT_EXTERNAL_DB_MAP",
    "MapFormatError",
    "get_external_db_map",
]

from importlib.resources import as_file, files
from pathlib import Path
from typing import Dict

# Provide the default map from a data file
default_map_res = files("ensembl.io.genomio.data.external_db_map").joinpath("default.txt")
with as_file(default_map_res) as default_map_path:
    DEFAULT_EXTERNAL_DB_MAP = default_map_path


class MapFormatError(ValueError):
    """Error when parsing the db map file."""


def get_external_db_map(map_file: Path) -> Dict:
    """Get an external_db map from a file.

    Returns:
        Dict with keys as alternate names, and values as standard name.

    """
    db_map: Dict[str, str] = {}
    with map_file.open("r") as map_fh:
        for line in map_fh:
            line = line.rstrip()
            if line.startswith("#") or line.startswith(" ") or line == "":
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise MapFormatError(f"External db file is not formatted correctly for: {line}")
            (main_name, alt_name) = parts[0:2]
            db_map[alt_name] = main_name
    return db_map
