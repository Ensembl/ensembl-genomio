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

import argschema
import jsonschema


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


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    json_file = argschema.fields.InputFile(required=True, metadata={"description": "JSON file to check"})
    json_schema = argschema.fields.InputFile(
        required=True, metadata={"description": "JSON schema to validate against"}
    )


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    json_schema_validator(mod.args["json_file"], mod.args["json_schema"])
