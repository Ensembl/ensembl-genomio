#!/usr/bin/env python3
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

import argschema

from ensembl.io.genomio.events_loader import EventCollection


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


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    # Server parameters
    input_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Input file from gene_diff"}
    )
    deletes_file = argschema.fields.InputFile(
        required=True,
        metadata={"description": "Deleted genes files (apart from the deletes from the gene diff)."},
    )
    output_file = argschema.fields.OutputFile(required=True, metadata={"description": "Formatted event file"})
    map_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Mapping tab file with 2 columns: old_id, new_id"}
    )
    release_name = argschema.fields.String(
        required=True, metadata={"description": "Name of the release for all event."}
    )
    release_date = argschema.fields.String(
        required=True, metadata={"description": "Date of the release for all events."}
    )


def main() -> None:
    """Main entrypoint"""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    args = mod.args

    # Start
    events = EventCollection()
    if args.get("deletes_file"):
        deleted_genes = load_list(args.get("deletes_file"))
        events.add_deletes(deleted_genes, args.get("release_name"), args.get("release_date"))
    events.load_events_from_gene_diff(
        args.get("input_file"), args.get("release_name"), args.get("release_date")
    )
    mapper = IdsMapper(args["map_file"])
    events.remap_to_ids(mapper.map)
    events.write_events_to_file(args["output_file"])


if __name__ == "__main__":
    main()
