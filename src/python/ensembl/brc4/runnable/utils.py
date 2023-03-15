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


import json
from pathlib import Path
from typing import Any


def print_json(path: Path, data: Any) -> None:
    """Generic data json dumper to a file.

    Args:
        path: Path to the json to create.
        data: Any data to store.
    """
    with path.open("w") as json_out:
        json_out.write(json.dumps(data, sort_keys=True, indent=4))


def get_json(json_path: Path) -> Any:
    """Generic data json loader.

    Args:
        path: Path to the json file to load.
    """
    with json_path.open("r") as json_file:
        return json.load(json_file)
