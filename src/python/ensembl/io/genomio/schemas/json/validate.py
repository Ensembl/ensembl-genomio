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

Examples:

    >>> from ensembl.io.genomio.schemas import json
    >>> json.schema_validator(json_file="functional_annotation.json", json_schema="functional_annotation")
    >>> json.schema_validator(json_file="functional_annotation.json", json_schema="genome")
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "ensembl-genomio/src/python/ensembl/io/genomio/schemas/json/validate.py", line 63,
      in schema_validator
        jsonschema.validate(instance=content, schema=schema)
      File ".venv/dev/lib/python3.10/site-packages/jsonschema/validators.py", line 1306, in validate
        raise error
    <list of all the elements from functional_annotation.json that failed validation>

"""

__all__ = ["schema_validator"]

from importlib.resources import as_file, files
import json
from os import PathLike
from pathlib import Path
from typing import Union

import jsonschema

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser


_JSON_SCHEMAS = {}
for schema_file in files("ensembl.io.genomio.data.schemas").iterdir():
    with as_file(schema_file) as file_path:
        if file_path.suffix == ".json":
            _JSON_SCHEMAS[file_path.stem] = file_path


def schema_validator(json_file: PathLike, json_schema: Union[str, PathLike]) -> None:
    """Validates a JSON file with the provided JSON schema.

    Args:
        json_file: Path to the JSON file to check.
        json_schema: JSON schema to validate `json_file` against, either a string matching a existing
            schema (in data/schemas) or a JSON schema file.

    """
    # Open IO for JSON files and validate it
    with Path(json_file).open("r") as fh:
        content = json.load(fh)
    # Find the json_schema file if a known identifier is provided (if not, treat it as a file path)
    if isinstance(json_schema, str) and (json_schema in _JSON_SCHEMAS):
        json_schema = _JSON_SCHEMAS[json_schema]
    with Path(json_schema).open("r") as fh:
        schema = json.load(fh)
    jsonschema.validate(instance=content, schema=schema)


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description="Validates a JSON file against a JSON schema.")
    parser.add_argument_src_path("--json_file", required=True, help="JSON file to check")
    parser.add_argument(
        "--json_schema", required=True, choices=_JSON_SCHEMAS.keys(), help="JSON schema to validate against"
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    args = parser.parse_args()

    schema_validator(args.json_file, args.json_schema)
