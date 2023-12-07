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
"""Extract gpad lines for a given list of IDs and species."""

__all__ = []

import logging
from pathlib import Path
from typing import Any, Dict, List
from os import PathLike

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args
from ensembl.io.genomio.utils import open_gz_file


_DATABASE = "UniProtKB"


def scan_gpad(gpad_file: PathLike, map_file: PathLike, output_file: PathLike) -> None:
    with open_gz_file(gpad_file) as gpad_fh, Path(output_file).open("w") as output_fh:
        logging.info("Get header")
        line_one = gpad_fh.readline()
        if line_one.startswith("gpa-version:"):
            output_fh.write(line_one)
        else:
            raise Exception(f"GPAD file must start with a header gpa-version. First line is {line_one}")

        logging.info("Get comments")
        for line in gpad_fh:
            if not line.startswith("!"):
                break
            if not (line.startswith("!Generated:") or line.startswith("!GO-version:")):
                continue
            output_fh.write(line)

        logging.info("Get entries")
        max_entries = 3
        n_entries = 0
        for line in gpad_fh:
            n_entries += 1
            logging.debug(f"Reading entry {n_entries}")
            if n_entries > max_entries:
                break

            parts = line.split("\t")
            if parts[0] != _DATABASE:
                logging.debug(f"Wrong db '{parts[0]}'")
                continue
            output_fh.write(line)


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description="Extract gpad entries from a list of IDs.")
    parser.add_argument_src_path("--map")
    parser.add_argument_src_path("--gpad")
    parser.add_argument_dst_path("--output")
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    scan_gpad(args.gpad, args.map, args.output)
