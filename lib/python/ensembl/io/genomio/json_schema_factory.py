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
"""Generates one JSON file per metadata type inside `manifest`, including the manifest itself.

Can be imported as a module and called as a script as well, with the same parameters and expected outcome.
"""

import ast
import json
import shutil
from pathlib import Path
from typing import List, Union

import argschema


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    manifest_dir = argschema.fields.InputDir(
        required=True, description="Folder containing the 'manifest.json' file to check"
    )
    metadata_types = argschema.fields.String(
        required=True, description="Metadata types to extract (in a list-like string)"
    )


def json_schema_factory(manifest_dir: Union[str, Path], metadata_types: List[str]) -> None:
    """Generates one JSON file per metadata type inside `manifest`, including "manifest.json" itself.

    Each JSON file will have the file name of the metadata type, e.g. "seq_region.json".

    Args:
        manifest_dir: Path to the folder with the manifest JSON file to check.
        metadata_types: Metadata types to extract from `manifest` as JSON files.

    """
    manifest_path = Path(manifest_dir, "manifest.json")
    with manifest_path.open() as manifest_file:
        content = json.load(manifest_file)
        shutil.copyfile(manifest_path, "manifest.json")
        # Use dir name from the manifest
        for name in content:
            if "file" in content[name]:
                file_name = content[name]["file"]
                content[name] = str(manifest_path.parent / file_name)
            else:
                for key in content[name]:
                    if "file" in content[name][key]:
                        file_name = content[name][key]["file"]
                        content[name][key] = str(manifest_path.parent / file_name)
        # Check the other JSON schemas
        for metadata_key in metadata_types:
            if metadata_key in content:
                shutil.copyfile(content[metadata_key], f"{metadata_key}.json")


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    # mod.args["metadata_types"] will be a list-like string that needs to be parsed to List[str]
    metadata_types = ast.literal_eval(mod.args["metadata_types"])
    json_schema_factory(mod.args["manifest_dir"], metadata_types)


# if __name__ == "__main__":
#     main()
