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


def compare_assembly(ncbi_main: Dict, ncbi_organella: Dict, core: Dict) -> Dict:

    # First count the organella
    ncbi_num_organella = len(ncbi_organella)
    core_num_organella = 0
    core_num_chrs = 0
    for loc, loc_count in core["locations"].items():
        if loc == "nuclear_chromosome":
            core_num_chrs += loc_count
        else:
            core_num_organella += loc_count
    
    print(f"Core Chromosomes: {core_num_chrs}")
    print(f"Core Organella: {core_num_organella}")

    # Our core stats count Organella chromosomes, sanity check here
    core_adjusted_chrs = core["coord_system"].get("chromosome", 0) - core_num_organella
    if core_adjusted_chrs != core_num_chrs:
        raise ValueError(f"Core stats chromosomes number doesn't add up: {core_adjusted_chrs} vs {core_num_chrs}")
    
    core_num_scaffolds = core["coord_system"].get("scaffold", 0)
    print(f"Core scaffolds: {core_num_scaffolds}")

    # NCBI includes the chromodomes in its stats
    core_adjusted_scaffolds = core_num_scaffolds + core_num_chrs
    print(f"Core adjusted scaffolds: {core_adjusted_scaffolds}")

    # Compare:
    # The number of scaffolds
    comp = {
        "scaffolds": _diff_count(ncbi_main["number_of_scaffolds"], core_adjusted_scaffolds),
        "chromosomes": _diff_count(ncbi_main["total_number_of_chromosomes"], core_adjusted_chrs),
        "organella": _diff_count(ncbi_num_organella, core_num_organella),
    }

    return comp


def _diff_count(ncbi_count: int, core_count: int) -> Dict:
    comp = {
        "ncbi": ncbi_count,
        "core": core_count,
        "diff": ncbi_count - core_count
    }
    return comp


def compare_annotation(ncbi: Dict, core: Dict) -> Dict:
    comp = {}
    return comp


def compare_stats(ncbi: Dict, core: Dict) -> Dict:
    """Compare stats from NCBI and our core."""

    ncbi_assembly_stats = ncbi.get("assembly_stats")
    ncbi_organella = ncbi.get("organelle_info")
    ncbi_annotation_stats = ncbi.get("annotation_stats")
    core_assembly_stats = core.get("assembly_stats")
    core_annotation_stats = core.get("annotation_stats")

    comp = {
        "assembly_diff": compare_assembly(ncbi_assembly_stats, ncbi_organella, core_assembly_stats),
    }
    if core_annotation_stats is not None and ncbi_annotation_stats is not None:
        comp["annotation_diff"] = compare_annotation(ncbi_annotation_stats, core_annotation_stats)

    new_core = {
        "core_stats": core,
        "ncbi_stats": {"assembly_stats": ncbi_assembly_stats, "anotation_stats": core_annotation_stats},
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
        ncbi_stats = json.load(ncbi_fh)["reports"][0]
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
