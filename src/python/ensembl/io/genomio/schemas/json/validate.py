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
"""Validates a JSON file with the provided JSON schema."""

__all__ = ["json_schema_validator"]

import json
from os import PathLike
from pathlib import Path

import jsonschema

from ensembl.utils.argparse import ArgumentParser


def json_schema_validator(json_file: PathLike, json_schema: PathLike) -> None:
    """Validates a JSON file with the provided JSON schema.

    Args:
        json_file: Path to the JSON file to check.
        json_schema: JSON schema to validate `json_file` against.

    Example:
        check_json_schema --json_file <JSON_FILE> --json_schema <JSON_SCHEMA>
    """

    # Open IO for JSON files and validate it
    with Path(json_file).open("r") as fh:
        content = json.load(fh)
    with Path(json_schema).open("r") as fh:
        schema = json.load(fh)
    jsonschema.validate(instance=content, schema=schema)


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description="Validates a JSON file against a JSON schema")
    parser.add_argument_src_path("--json_file", required=True, help="JSON file to check")
    parser.add_argument_src_path("--json_schema", required=True, help="JSON schema to validate against")
    args = parser.parse_args()

    json_schema_validator(**vars(args))
