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
"""Creates a manifest file in a folder depending on the file names ends."""

__all__ = ["Manifest"]

import hashlib
import logging
import json
from pathlib import Path

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


class MannifestError(Exception):
    """Could not load a manifest file."""


class Manifest:
    """Records of a manifest file and its files and md5 checksums."""

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
    alias_names = {
        "gene_models": "gff3",
        "dna": "fasta_dna",
        "pep": "fasta_pep",
    }
    names = {name: name for name in same_names}
    names = {**names, **alias_names}
    multi_files = {"agp"}

    def __init__(self, manifest_dir: Path) -> None:
        """Initializes a manifest with the directory containing the files (and a manifest if it exists).
        
        Args:
            manifest_dir: directory where the files are contained.
        """
        self.dir = manifest_dir
        self.path = manifest_dir / "manifest.json"
        self.files = {}

    def create(self) -> Path:
        """Creates a manifest file from the files in a directory."""
        self.get_files_checksums()
        with self.path.open("w") as json_out:
            json_out.write(json.dumps(self.files, sort_keys=True, indent=4))

    def get_files_checksums(self) -> None:
        """Records all the files in the directory with their checksum."""
        manifest_files = {}
        for subfile in self.dir.iterdir():
            logging.debug(f"Check file {subfile} ({subfile.stem}, {subfile.suffix})")
            used_file = False
            if subfile.is_dir():
                logging.warning("Can't create manifest for subdirectory")
                continue

            # Delete and skip empty files
            if subfile.stat().st_size == 0:
                logging.warning(f"Skip and delete empty file: {subfile}")
                subfile.unlink()
                continue

            for name, standard_name in self.names.items():
                # Either the last element of the stem or the suffix is a known name
                if subfile.stem.endswith(name) or subfile.suffix == f".{name}":
                    logging.debug(f"Matched to {name} ({standard_name}) = {subfile}")
                    used_file = True
                    md5 = self._get_md5sum(subfile)
                    file_obj = {"file": subfile.name, "md5sum": md5}

                    # Multiple files stored, each with a name
                    if standard_name in self.multi_files:
                        manifest_files.setdefault(standard_name, {})
                        obj_name = self._prepare_object_name(subfile, name, manifest_files[standard_name])
                        manifest_files[standard_name][obj_name] = file_obj

                    # Single file/init
                    else:
                        manifest_files[standard_name] = file_obj

            if not used_file:
                logging.warning(f"File {subfile} was not included in the manifest")

        self.files = manifest_files

    def _prepare_object_name(self, subfile: Path, name: str, manifest_dict: dict) -> str:
        # Prepare object name
        obj_name = "file"
        try:
            # If we recognize the suffix, then the name is the part after the last "_"
            if subfile.suffix == f".{name}":
                obj_name = subfile.stem.split(sep="_")[-1]
            # If we recognize the end of the name, then the name is the part before the last "_"
            else:
                obj_name = subfile.stem.split(sep="_")[-2]
        except IndexError:
            pass

        # Add number if duplicate name
        obj_name_base = obj_name
        count = 1
        while obj_name in manifest_dict:
            obj_name = f"{obj_name_base}.{count}"
            count += 1
            if count >= 10:
                raise ValueError(f"Too many files with same name {obj_name_base}")
        return obj_name

    def load(self) -> None:
        """Load the content of an existing manifest file."""
        if not self.path.exists():
            raise MannifestError(f"Cannot load non-existing manifest file: {self.path}")

        with self.path.open("r") as manifest_fh:
            manifest = json.load(manifest_fh)

            # Use dir name from the manifest
            for name in manifest:
                if "file" in manifest[name]:
                    file_path = self.dir / manifest[name]["file"]
                    # check if the md5sum is correct
                    md5sum = manifest[name]["md5sum"]
                    self._check_md5sum(file_path, md5sum)

                    manifest[name] = file_path
                else:
                    for f in manifest[name]:
                        if "file" in manifest[name][f]:
                            file_path = self.dir / manifest[name][f]["file"]
                            # check if the md5sum is correct
                            md5sum = manifest[name][f]["md5sum"]
                            self._check_md5sum(file_path, md5sum)

                            manifest[name][f] = file_path

    @staticmethod
    def _get_md5sum(file_path: Path) -> str:
        """Returns the md5 checksum for a given file."""
        with file_path.open("rb") as f:
            data_bytes = f.read()
            return hashlib.md5(data_bytes).hexdigest()

    def _check_md5sum(self, file_path: Path, md5sum: str) -> None:
        """Checks a file against an md5 checksum, raises a ManifestError if the checksum fails.

        Args:
            file_path: Path to a genome file.
            md5sum: MD5 hash for the files.
        """
        if self._get_md5sum(file_path) != md5sum:
            raise MannifestError(f"Invalid md5 checksum for {file_path}")


def main() -> None:
    """Main entrypoint."""
    parser = ArgumentParser(
        description="Compare the genomic data between the files present in a manifest file."
    )
    parser.add_argument_dst_path(
        "--manifest_dir", required=True, help="Folder where to create a manifest file"
    )
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    manifest = Manifest(args.manifest_dir)
    manifest.create()


if __name__ == "__main__":
    main()
