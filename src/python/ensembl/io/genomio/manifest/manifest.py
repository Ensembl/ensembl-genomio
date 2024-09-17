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
"""Representation of a manifest file."""

__all__ = ["Manifest", "ManifestError"]

import hashlib
import json
import logging
from pathlib import Path
from typing import Any, TypeAlias


ManifestDict: TypeAlias = dict[str, dict[str, Any]]


class ManifestError(Exception):
    """Could not load a manifest file."""


class Manifest:
    """Records of a manifest file and its files and md5 checksums."""

    _same_names = {
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
    _alias_names = {
        "gene_models": "gff3",
        "dna": "fasta_dna",
        "pep": "fasta_pep",
    }
    _same_names_dict = {name: name for name in _same_names}
    names = {**_same_names_dict, **_alias_names}
    multi_files = {"agp"}

    def __init__(self, manifest_dir: Path) -> None:
        """Initializes a manifest with the directory containing the files (and a manifest if it exists).

        Args:
            manifest_dir: directory where the files are contained.
        """
        self.root_dir = manifest_dir
        self.file_path = manifest_dir / "manifest.json"
        self.files: dict = {}

    def create(self) -> None:
        """Creates a manifest file from the files in a directory."""
        self.get_files_checksums()
        with self.file_path.open("w") as json_out:
            json_out.write(json.dumps(self.files, sort_keys=True, indent=4))

    def get_files_checksums(self) -> ManifestDict:
        """Records all the files in the directory with their checksum."""
        manifest_files: ManifestDict = {}
        for subfile in self.root_dir.iterdir():
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
        return self.files

    def _prepare_object_name(
        self, subfile: Path, name: str, manifest_file_dict: dict[str, dict[str, str]]
    ) -> str:
        # Prepare object name
        try:
            # If we recognize the suffix, then the name is the part after the last "_"
            if subfile.suffix == f".{name}":
                obj_name = subfile.stem.split(sep="_")[-1]
            # If we recognize the end of the name, then the name is the part before the last "_"
            else:
                obj_name = subfile.stem.split(sep="_")[-2]
        except IndexError:
            obj_name = "file"

        # Add number if duplicate name
        obj_name_base = obj_name
        count = 1
        while obj_name in manifest_file_dict.keys():
            obj_name = f"{obj_name_base}.{count}"
            count += 1
            if count >= 10:
                raise ValueError(f"Too many files with same name {obj_name_base}")
        return obj_name

    def load(self) -> ManifestDict:
        """Load the content of an existing manifest file."""
        if not self.file_path.exists():
            raise ManifestError(f"Cannot load non-existing manifest file: {self.file_path}")

        with self.file_path.open("r") as manifest_fh:
            manifest = json.load(manifest_fh)

            # Use dir name from the manifest
            for name in manifest:
                if "file" in manifest[name]:
                    file_path = self.root_dir / manifest[name]["file"]
                    # check if the md5sum is correct
                    md5sum = manifest[name]["md5sum"]
                    self._check_md5sum(file_path, md5sum)
                else:
                    for f in manifest[name]:
                        file_path = self.root_dir / manifest[name][f]["file"]
                        # check if the md5sum is correct
                        md5sum = manifest[name][f]["md5sum"]
                        self._check_md5sum(file_path, md5sum)

            self.files = manifest
        return self.files

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
        file_md5sum = self._get_md5sum(file_path)
        if file_md5sum != md5sum:
            raise ManifestError(f"Invalid md5 checksum for {file_path}: got {file_md5sum}, expected {md5sum}")
