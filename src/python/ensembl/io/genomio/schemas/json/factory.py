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

__all__ = ["schema_factory"]

import json
from os import PathLike
from pathlib import Path
import shutil

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def schema_factory(manifest_dir: PathLike) -> None:
    """Generates one JSON file per metadata type inside `manifest`, including "manifest.json" itself.

    Each JSON file will have the file name of the metadata type, e.g. "seq_region.json".

    Args:
        manifest_dir: Path to the folder with the manifest JSON file to check.
    """
    metadata_types = ["seq_region", "genome", "functional_annotation"]
    manifest_path = Path(manifest_dir)
    with manifest_path.open() as manifest_file:
        content = json.load(manifest_file)
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
                        validate_files= [filepath, metadata_key]
                else:
                    validate_files = [json_files[metadata_key], metadata_key]
        return validate_files

def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Generates one JSON file per metadata type in the provided manifest, including itself."
    )
    parser.add_argument_src_path(
        "--manifest_dir", required=True, help="Folder containing the 'manifest.json' file to check"
    )
    args = parser.parse_args()
    init_logging_with_args(args)

    schema_factory(**vars(args))
