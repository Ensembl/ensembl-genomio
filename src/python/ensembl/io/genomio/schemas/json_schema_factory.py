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
"""Generates one JSON file per metadata type inside `manifest`, including the manifest itself."""

__all__ = ["json_schema_factory"]

import ast
import json
from os import PathLike
from pathlib import Path
import shutil
from typing import List

import argschema


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    manifest_dir = argschema.fields.InputDir(
        required=True, metadata={"description": "Folder containing the 'manifest.json' file to check"}
    )
    metadata_types = argschema.fields.String(
        required=True, metadata={"description": "Metadata types to extract (in a list-like string)"}
    )
    output_dir = argschema.fields.InputDir(
        required=False, dump_default=".", metadata={"description": "Folder to store the produced files"}
    )


def json_schema_factory(manifest_dir: PathLike, metadata_types: List[str], output_dir: PathLike) -> None:
    """Generates one JSON file per metadata type inside `manifest`, including "manifest.json" itself.

    Each JSON file will have the file name of the metadata type, e.g. "seq_region.json".

    Args:
        manifest_dir: Path to the folder with the manifest JSON file to check.
        metadata_types: Metadata types to extract from `manifest` as JSON files.
        output_dir: Path to the folder where to generate the JSON files.

    """
    manifest_path = Path(manifest_dir, "manifest.json")
    with manifest_path.open() as manifest_file:
        content = json.load(manifest_file)
        shutil.copyfile(manifest_path, Path(output_dir, "manifest.json"))
        json_files = {}
        # Use dir name from the manifest
        for name in content:
            if "file" in content[name]:
                file_name = content[name]["file"]
                json_files[name] = manifest_path.parent / file_name
            else:
                for key in content[name]:
                    if "file" in content[name][key]:
                        file_name = content[name][key]["file"]
                        json_files[name] = {key: manifest_path.parent / file_name}
        # Check the other JSON schemas
        for metadata_key in metadata_types:
            if metadata_key in json_files:
                if isinstance(json_files[metadata_key], dict):
                    for key, filepath in json_files[metadata_key].items():
                        shutil.copyfile(filepath, Path(output_dir, f"{metadata_key}_{key}.json"))
                else:
                    shutil.copyfile(json_files[metadata_key], Path(output_dir, f"{metadata_key}.json"))


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    # mod.args["metadata_types"] will be a list-like string that needs to be parsed to List[str]
    metadata_types = ast.literal_eval(mod.args["metadata_types"])
    json_schema_factory(mod.args["manifest_dir"], metadata_types, mod.args["output_dir"])
