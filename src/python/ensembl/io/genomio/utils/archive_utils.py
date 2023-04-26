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
"""TODO"""

__all__ = ["SUPPORTED_ARCHIVE_FORMATS", "extract_file"]

import gzip
from os import PathLike
from pathlib import Path
import shutil

import argschema


# Each registered format is a tuple, `(name, extensions, description)`
SUPPORTED_ARCHIVE_FORMATS = [ext for elem in shutil.get_unpack_formats() for ext in elem[1]]


def extract_file(src_file: PathLike, dst_dir: PathLike) -> None:
    """Extracts the `src_file` into `dst_dir`.

    If the file is not an archive, it will only be copied to `dst_dir`.

    Args:
        src_file: Path to the file to unpack.
        dst_dir: Path to where extract the file.

    """
    extension = Path(src_file).suffix
    file_base = src_file.strip(extension)
    final_path = Path(dst_dir) / file_base

    if extension == ".gz":
        with gzip.open(src_file, "rb") as f_in:
            with open(final_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    elif extension in SUPPORTED_ARCHIVE_FORMATS:
        shutil.unpack_archive(src_file, dst_dir)
    else:
        shutil.copy(src_file, dst_dir)


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by the entry point of `extract_file_cli`."""

    src_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Path to the file to unpack"}
    )
    dst_dir = argschema.fields.OutputDir(
        required=True, 
        dump_default=".", 
        metadata={"description": "Output folder to where extract the file"}
    )


def extract_file_cli() -> None:
    """Entry-point for the `extract_file` method."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    extract_file(mod.args["src_file"], mod.args["dst_dir"])
