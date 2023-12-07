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
_ECO_TO_GO = {
    "ECO:0000269": "EXP",
    "ECO:0007005": "HDA",
    "ECO:0007007": "HEP",
    "ECO:0007003": "HGI",
    "ECO:0007001": "HMP",
    "ECO:0006056": "HTP",
    "ECO:0000318": "IBA",
    "ECO:0000319": "IBD",
    "ECO:0000305": "IC",
    "ECO:0000314": "IDA",
    "ECO:0000501": "IEA",
    "ECO:0000256": "IEA",
    "ECO:0007322": "IEA",
    "ECO:0000265": "IEA",
    "ECO:0000249": "IEA",
    "ECO:0000363": "IEA",
    "ECO:0000270": "IEP",
    "ECO:0000317": "IGC",
    "ECO:0000354": "IGC",
    "ECO:0000316": "IGI",
    "ECO:0000320": "IKR",
    "ECO:0000315": "IMP",
    "ECO:0000353": "IPI",
    "ECO:0000321": "IRD",
    "ECO:0000247": "ISA",
    "ECO:0000255": "ISM",
    "ECO:0000266": "ISO",
    "ECO:0000250": "ISS",
    "ECO:0000031": "ISS",
    "ECO:0000303": "NAS",
    "ECO:0000307": "ND",
    "ECO:0000245": "RCA",
    "ECO:0000304": "TAS",
}


def get_mapping(map_file: PathLike) -> Dict[str, Any]:

    mapping = {}
    with open(map_file, "r") as map_fh:
        for line in map_fh:
            if line == "":
                continue
            (db_id, _, production_name) = line.strip().split("\t")

            if db_id in mapping:
                mapping[db_id].append(production_name)
            else:
                mapping[db_id] = [production_name]
    return mapping


def scan_gpad(gpad_file: PathLike, map_file: PathLike, output_dir: PathLike) -> None:
    mapping = get_mapping(map_file)

    with open_gz_file(gpad_file) as gpad_fh:
        header = []
        done = set()

        # Get the header of all files
        logging.debug("Get header")
        line_one = gpad_fh.readline()
        if line_one.startswith("gpa-version:"):
            header.append(line_one)
        else:
            raise ValueError(f"GPAD file must start with a header gpa-version. First line is {line_one}")

        logging.debug("Get comments")
        for line in gpad_fh:
            if not line.startswith("!"):
                break
            if not (line.startswith("!Generated:") or line.startswith("!GO-version:")):
                continue
            header.append(line)

        logging.debug("Get entries")
        for line in gpad_fh:

            parts = line.split("\t")
            db = parts[0]
            if db != _DATABASE:
                continue

            accession = parts[1]
            if accession not in mapping:
                continue

            eco = parts[5]
            go_evidence = ""
            if eco in _ECO_TO_GO:
                go_evidence = _ECO_TO_GO[eco]

            prod_names = mapping[accession]
            for production_name in prod_names:
                parts[-1] = f"tgt_species={production_name}|go_evidence={go_evidence}\n"
                with open(output_dir / f"{production_name}.gpad", "a") as output_fh:
                    if production_name not in done:
                        for header_line in header:
                            output_fh.write(header_line)
                        done.add(production_name)
                    output_fh.write("\t".join(parts))


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description="Extract gpad entries from a list of IDs.")
    parser.add_argument_src_path("--map")
    parser.add_argument_src_path("--gpad")
    parser.add_argument_dst_path("--output_dir")
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    scan_gpad(args.gpad, args.map, args.output_dir)
