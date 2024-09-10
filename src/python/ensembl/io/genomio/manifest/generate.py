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

__all__ = ["ManifestMaker"]

import hashlib
import logging
import json
from pathlib import Path

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


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
    alias_names = {
        "gene_models": "gff3",
        "dna": "fasta_dna",
        "pep": "fasta_pep",
    }
    names = {name: name for name in same_names}
    names = {**names, **alias_names}
    multi_files = {"agp"}

    def __init__(self, manifest_dir: Path) -> None:
        self.dir = manifest_dir

    def create_manifest(self) -> Path:
        """Create the manifest file."""
        manifest_data = self.get_files_checksums()
        manifest_path = self.dir / "manifest.json"
        with manifest_path.open("w") as json_out:
            json_out.write(json.dumps(manifest_data, sort_keys=True, indent=4))
        return manifest_path

    def get_files_checksums(self) -> dict[str, dict]:
        """Compute the checksum of all the files in the directory."""
        manifest_files: dict[str, dict] = {}
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

        return manifest_files

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

    @staticmethod
    def _get_md5sum(file_path: Path) -> str:
        with file_path.open("rb") as f:
            data_bytes = f.read()
            return hashlib.md5(data_bytes).hexdigest()


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

    maker = ManifestMaker(args.manifest_dir)
    maker.create_manifest()


if __name__ == "__main__":
    main()
