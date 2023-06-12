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
"""Creates a manifest file in a folder depending on the file names ends.
"""

import hashlib
import json
from pathlib import Path

import argschema


class ManifestMaker:
    """Given a directory with genomic files, create a manifest json file for them."""

    same_names = {
        "gff3",
        "fasta_dna",
        "fasta_pep",
        "functional_annotation",
        "genome",
        "seq_attrib",
        "seq_region",
        "agp",
        "events",
    }
    alias_names = {"gene_models": "gff3", "dna": "fasta_dna", "pep": "fasta_pep"}
    names = {name: name for name in same_names}
    names = {**names, **alias_names}

    def __init__(self, manifest_dir: Path) -> None:
        self.dir = manifest_dir

    def create_manifest(self) -> Path:
        """Create the manifest file."""
        manifest_data = self.get_files_checksums()
        manifest_path = self.dir / "manifest.json"
        with manifest_path.open("w") as json_out:
            json_out.write(json.dumps(manifest_data, sort_keys=True, indent=4))
        return manifest_path

    def get_files_checksums(self):
        """Compute the checksum of all the files in the directory."""
        manifest_files = {}
        for subfile in self.dir.iterdir():
            used_file = False
            if subfile.is_dir():
                print("Can't create manifest for subdirectory")
                continue

            for name, standard_name in self.names.items():
                if subfile.stem.endswith(name):
                    used_file = True
                    md5 = self._get_md5sum(subfile)
                    file_obj = {"file": subfile.name, "md5sum": md5}
                    if standard_name in manifest_files:
                        if isinstance(manifest_files[standard_name], list):
                            manifest_files[standard_name].append(file_obj)
                        else:
                            # Convert to a list
                            manifest_files[standard_name] = [manifest_files[standard_name], file_obj]

                    else:
                        manifest_files[standard_name] = file_obj
                    break

            if not used_file:
                print(f"File {subfile} was not included in the manifest")

        return manifest_files

    @staticmethod
    def _get_md5sum(file_path: Path) -> str:
        with file_path.open("rb") as f:
            data_bytes = f.read()
            return hashlib.md5(data_bytes).hexdigest()


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    # Server parameters
    manifest_dir = argschema.fields.files.InputDir(
        required=True, metadata={"description": "Folder where to create a manifest file."}
    )


def main() -> None:
    """Main entrypoint."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)

    maker = ManifestMaker(Path(mod.args["manifest_dir"]))
    manifest_path = maker.create_manifest()

    if mod.args.get("output_json"):
        out_data = {"manifest_path": str(manifest_path)}
        with Path(mod.args["output_json"]).open("w") as output_fh:
            output_fh.write(json.dumps(out_data))


if __name__ == "__main__":
    main()
