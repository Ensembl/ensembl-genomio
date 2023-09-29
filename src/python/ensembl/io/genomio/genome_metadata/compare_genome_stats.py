#!/usr/bin/env python
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
"""Compare stats in a JSON from NCBI dataset and a json from our core.
Returns the json from our core including a section with comparisons.
"""

import json
from pathlib import Path
from typing import Any, Dict

import argschema


def compare_assembly(ncbi: Dict, core: Dict) -> Dict:
    comp = {}
    return comp


def compare_annotation(ncbi: Dict, core: Dict) -> Dict:
    comp = {}
    return comp


def compare_stats(ncbi: Dict, core: Dict) -> Dict:
    """Compare stats from NCBI and our core."""
    comp = {
        "assembly_diff": compare_assembly(ncbi, core),
        "annotation_diff": compare_annotation(ncbi, core)
    }

    new_core = {
        "core_stats": core,
        "ncbi_comparison": comp
    }
    return new_core


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    # Server parameters
    ncbi_stats = argschema.fields.InputFile(required=True, metadata={"description": "NCBI json file"})
    core_stats = argschema.fields.InputFile(required=True, metadata={"description": "Core json file"})


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    args = mod.args

    with open(args["ncbi_stats"]) as ncbi_fh:
        ncbi_stats = json.load(ncbi_fh)
    with open(args["core_stats"]) as core_fh:
        core_stats = json.load(core_fh)
    all_stats = compare_stats(ncbi_stats, core_stats)

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(all_stats, indent=2, sort_keys=True))
    else:
        print(json.dumps(all_stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
