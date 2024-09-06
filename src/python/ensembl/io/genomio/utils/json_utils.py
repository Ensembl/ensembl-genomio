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
"""Utils to deal with JSON files."""

__all__ = ["get_json", "print_json"]

import json
from pathlib import Path
from typing import Any

from ensembl.utils import StrPath


def get_json(src_path: StrPath, **kwargs: Any) -> Any:
    """Generic data JSON loader.

    Args:
        src_path: Path to the JSON file to load.

    """
    with Path(src_path).open("r", encoding="utf-8") as json_file:
        return json.load(json_file, **kwargs)


def print_json(dst_path: StrPath, data: Any, **kwargs: Any) -> None:
    """Generic data JSON dumper to a file, with keys sorted and pretty-printed with indent 4 by default.

    Args:
        dst_path: Path to the JSON file to create.
        data: Any data to store into the file.

    """
    kwargs.setdefault("sort_keys", True)
    kwargs.setdefault("indent", 4)
    with Path(dst_path).open("w", encoding="utf-8") as json_file:
        json_file.write(json.dumps(data, **kwargs))
