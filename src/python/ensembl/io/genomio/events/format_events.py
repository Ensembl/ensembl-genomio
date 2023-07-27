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


from os import PathLike
from pathlib import Path
from typing import Dict

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
                (from_id, to_id) = line.split("\t")[0:2]
                mapping[from_id] = to_id

        return mapping


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    # Server parameters
    input_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Input file from gene_diff"}
    )
    output_file = argschema.fields.OutputFile(required=True, metadata={"description": "Formatted event file"})
    map_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Mapping tab file with 2 columns: old_id, new_id"}
    )


def main() -> None:
    """Main entrypoint"""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    args = mod.args

    # Start
    events = EventCollection()
    events.load_events_from_gene_diff(args["input_file"])
    mapper = IdsMapper(args["map_file"])
    events.remap_to_ids(mapper.map)
    events.write_events_to_file(args["output_file"])


if __name__ == "__main__":
    main()
