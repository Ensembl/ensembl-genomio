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
"""Module to map stable ids in a file, given a mapping."""

__all__ = ["IdsMapper", "load_list"]

from os import PathLike
from pathlib import Path
import re
from typing import Dict, List

import ensembl.io.genomio
from ensembl.io.genomio.events.load import EventCollection
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


class IdsMapper:
    """Simple mapper object, to cleanly get a mapping dict."""

    def __init__(self, map_file: PathLike) -> None:
        self.map = self._load_mapping(Path(map_file))

    def _load_mapping(self, map_file: Path) -> Dict[str, str]:
        """Return a mapping in a simple dict from a tab file with 2 columns: from_id, to_id.

        Args:
            map_file: Tab file path.
        """
        mapping = {}
        with map_file.open("r") as map_fh:
            for line in map_fh:
                if line == "":
                    continue
                items = line.split("\t")
                if len(items) < 2:
                    raise ValueError(f"Not 2 elements in {line}")
                (from_id, to_id) = items[0:2]
                mapping[from_id] = to_id

        return mapping


def load_list(list_file: Path) -> List[str]:
    """Return a simple list from a file."""
    items = set()
    empty_spaces = re.compile(r"\s+")
    with Path(list_file).open("r") as map_fh:
        for line in map_fh:
            line = re.sub(empty_spaces, "", line)
            if line == "":
                continue
            items.add(line)

    return list(items)


def main() -> None:
    """Main entrypoint"""
    parser = ArgumentParser(description="Map stable IDs in a file and produce an events file.")
    parser.add_argument_src_path("--input_file", required=True, help="Input file from gene_diff")
    parser.add_argument_src_path(
        "--deletes_file", required=True, help="Deleted genes file (apart from the deletes from the gene diff)"
    )
    parser.add_argument_src_path(
        "--map_file", required=True, help="Mapping tab file with 2 columns: old_id, new_id"
    )
    parser.add_argument("--release_name", required=True, metavar="NAME", help="Release name for all events")
    parser.add_argument("--release_date", required=True, metavar="DATE", help="Release date for all events")
    parser.add_argument_dst_path("--output_file", required=True, help="Output formatted event file")
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    events = EventCollection()
    deleted_genes = load_list(args.deletes_file)
    events.add_deletes(deleted_genes, args.release_name, args.release_date)
    events.load_events_from_gene_diff(args.input_file, args.release_name, args.release_date)
    mapper = IdsMapper(args.map_file)
    events.remap_to_ids(mapper.map)
    events.write_events_to_file(args.output_file)


if __name__ == "__main__":
    main()
