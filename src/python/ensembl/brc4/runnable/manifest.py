#!env python3
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


import hashlib
from pathlib import Path
import shutil
from typing import Dict

import eHive

from ensembl.io.genomio.utils.json_utils import print_json, get_json


class manifest(eHive.BaseRunnable):
    def run(self):
        manifest_files = self.param_required("manifest_files")
        output_dir = Path(self.param_required("output_dir"))

        # Assume there is a genome file and get metadata from it
        if not "genome" in manifest_files:
            genome_file = self.param("genome_json")
            if genome_file:
                manifest_files["genome"] = genome_file
            else:
                raise Exception("Processed genome metadata file is missing")

        # to configure the file names
        genome_data = self.get_genome_data(manifest_files)
        print(genome_data)

        # Get BRC4 component, organism if any
        component = None
        genome_name = self.param("genome_name")
        species = self.param("species")

        if not genome_name:
            if "BRC4" in genome_data:
                if "component" in genome_data["BRC4"]:
                    component = genome_data["BRC4"]["component"]
                if "organism_abbrev" in genome_data["BRC4"]:
                    genome_name = genome_data["BRC4"]["organism_abbrev"]

            # RETROCOMPATIBILITY
            if (
                genome_name == None
                and "species" in genome_data
                and "BRC4_organism_abbrev" in genome_data["species"]
            ):
                genome_name = genome_data["species"]["BRC4_organism_abbrev"]

            # Get the species name (production_name)
            if "species" in genome_data and "production_name" in genome_data["species"]:
                species = genome_data["species"]["production_name"]

            # Default species name
            if genome_name == None:
                if species != None:
                    genome_name = species
                else:
                    raise Exception("No species")

        # Define genome directory
        genome_dir = output_dir
        if component != None:
            genome_dir = genome_dir / component
        genome_dir = genome_dir / genome_name

        if not genome_dir.is_dir():
            genome_dir.mkdir(parents=True)
        print(f"Genome directory: '{genome_dir}'")

        # Move all files to the output dir, with a new name
        print(manifest_files)
        final_files = self.copy_files(genome_dir, manifest_files, species, genome_name)

        # Create th manifest file itself
        manifest_path = self.create_manifest(genome_dir, final_files)

        self.dataflow({"manifest": str(manifest_path)}, 2)

    def get_genome_data(self, files):
        if "genome" in files:
            return get_json(Path(files["genome"]))
        else:
            raise Exception("No genome file")

    def copy_files(self, genome_dir: Path, files: Dict, species: str, genome_name: str) -> Dict:
        final_files = {}
        for name, value in files.items():
            final_files[name] = {}

            # Nested files
            if isinstance(value, dict):
                for subname, subvalue in files[name].items():
                    file_path = Path(subvalue)
                    final_file = self.copy_file(file_path, genome_dir, species, genome_name)
                    if final_file:
                        final_files[name][subname] = final_file

            # Not nested files
            else:
                file_path = Path(value)
                final_file = self.copy_file(file_path, genome_dir, species, genome_name)
                if final_file:
                    final_files[name] = final_file

            if len(final_files[name]) == 0:
                del final_files[name]

        return final_files

    def copy_file(self, origin_path: Path, dir_path: Path, species: str, genome_name: str) -> Dict:
        if not origin_path:
            return None

        # Define the new file name
        file_name = origin_path.name
        if species != None and genome_name != None:
            file_name = file_name.replace(species, genome_name)
        final_path = dir_path / file_name
        print(f"{origin_path} -> {final_path} ({dir_path})")

        # Copy file
        shutil.copyfile(origin_path, final_path)

        # Compute md5sum
        md5sum = self.get_md5sum(final_path)

        file_data = {"file": file_name, "md5sum": md5sum}

        return file_data

    def get_md5sum(self, file_path: Path) -> str:
        with file_path.open("rb") as f:
            bytes = f.read()
            return hashlib.md5(bytes).hexdigest()

    def create_manifest(self, genome_dir: Path, files: Dict) -> Path:
        manifest_path = genome_dir / "manifest.json"
        print_json(manifest_path, files)
        return manifest_path
