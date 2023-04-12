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
"""TODO"""

__all__ = ["get_json", "print_json"]

import json
from os import PathLike
from pathlib import Path
import shutil
from typing import Any


def get_json(json_path: PathLike) -> Any:
    """Generic data JSON loader.

    Args:
        path: Path to the JSON file to load.

    """
    with Path(json_path).open("r") as json_file:
        return json.load(json_file)


def print_json(json_path: PathLike, data: Any) -> None:
    """Generic data JSON dumper to a file.

    Args:
        path: Path to the JSON to create.
        data: Any data to store into the file.

    """
    with Path(json_path).open("w") as json_file:
        json_file.write(json.dumps(data, sort_keys=True, indent=4))
