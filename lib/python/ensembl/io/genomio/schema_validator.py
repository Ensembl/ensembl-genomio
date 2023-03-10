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
"""Validates a JSON file with the provided JSON schema.

Can be imported as a module and called as a script as well, with the same parameters and expected outcome.
"""

import json
import jsonschema

import argschema


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    json_file = argschema.fields.InputFile(required=True, description="JSON file to check")
    json_schema = argschema.fields.InputFile(required=True, description="JSON schema to validate against")


def schema_validator(json_file: str, json_schema: str) -> None:
    """Validates a JSON file with the provided JSON schema.

    Args:
        json_file: Path to the JSON file to check.
        json_schema: JSON schema to validate `json_file` against.

    """
    with open(json_file) as fh:
        content = json.load(fh)
    with open(json_schema) as fh:
        schema = json.load(fh)
    jsonschema.validate(instance=content, schema=schema)


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    schema_validator(mod.args["json_file"], mod.args["json_schema"])


if __name__ == "__main__":
    main()
