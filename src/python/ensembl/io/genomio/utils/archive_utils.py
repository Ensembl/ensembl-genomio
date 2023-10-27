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
"""Utils to deal with archive files."""

__all__ = ["SUPPORTED_ARCHIVE_FORMATS", "extract_file"]

from contextlib import contextmanager
import gzip
from os import PathLike
from pathlib import Path
import shutil
from typing import Generator, TextIO

import argschema


# Each registered format is a tuple, `(name, extensions, description)`
SUPPORTED_ARCHIVE_FORMATS = [ext for elem in shutil.get_unpack_formats() for ext in elem[1]]


@contextmanager
def open_gz_file(file_path: PathLike) -> Generator[TextIO, None, None]:
    """Yields an open file object, even if the file is compressed with gzip.

    The file is expected to contain a text, and this can be used with the usual "with".

    Args:
        file_path: A file path to open.

    """
    this_file = Path(file_path)
    if this_file.suffix == ".gz":
        with gzip.open(this_file, "rt") as fh:
            yield fh
    else:
        with this_file.open("rt") as fh:
            yield fh


def extract_file(src_file: PathLike, dst_dir: PathLike) -> None:
    """Extracts the `src_file` into `dst_dir`.

    If the file is not an archive, it will be copied to `dst_dir`. `dst_dir` will be created if it
    does not exist.

    Args:
        src_file: Path to the file to unpack.
        dst_dir: Path to where extract the file.

    """
    src_file = Path(src_file)
    extensions = {"".join(src_file.suffixes[i:]) for i in range(0, len(src_file.suffixes))}
    file_base = src_file.stem
    dst_dir = Path(dst_dir)
    final_path = dst_dir / file_base

    if extensions.intersection(SUPPORTED_ARCHIVE_FORMATS):
        shutil.unpack_archive(src_file, dst_dir)
    else:
        # Replicate the functionality of shutil.unpack_archive() by creating `dst_dir`
        dst_dir.mkdir(parents=True, exist_ok=True)
        if extensions.intersection([".gz"]):
            with gzip.open(src_file, "rb") as f_in:
                with final_path.open("wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy(src_file, dst_dir)


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by the entry point of `extract_file_cli`."""

    src_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Path to the file to unpack"}
    )
    dst_dir = argschema.fields.OutputDir(
        required=True, dump_default=".", metadata={"description": "Output folder to where extract the file"}
    )


def extract_file_cli() -> None:
    """Entry-point for the `extract_file` method."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    extract_file(mod.args["src_file"], mod.args["dst_dir"])
