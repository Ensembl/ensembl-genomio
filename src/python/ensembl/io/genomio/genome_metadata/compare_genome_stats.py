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
import re
from typing import Any, Dict

import argschema


def _diff_dicts(ncbi: Dict[str, int], core: Dict[str, int]) -> Dict:
    diff = {}
    same = {}
    for key in ncbi.keys():
        if ncbi[key] == 0 and core[key] == 0:
            continue
        if ncbi[key] == core[key]:
            same[key] = ncbi[key]
            continue
        diff[key] = {"ncbi": ncbi[key], "core": core[key], "diff": core[key] - ncbi[key]}

    comp = {}
    if same:
        comp["same"] = same
    if diff:
        comp["different"] = diff
    return comp


def compare_assembly(ncbi: Dict, core: Dict) -> Dict:
    """Returns a compilation of count comparisons.
    Each comparison is a dict with the value from NCBI, from Core, and their diff.
    """

    # Prepare counts to be comparable to the NCBI stats
    ncbi_main = ncbi.get("assembly_stats", {})
    ncbi_info = ncbi.get("assembly_info", {})
    ncbi_organella = ncbi.get("organelle_info", [])

    # First count the organella
    core_num_organella = 0
    core_num_chrs = 0
    for loc, loc_count in core["locations"].items():
        if loc == "nuclear_chromosome":
            core_num_chrs += loc_count
        else:
            core_num_organella += loc_count

    # Our core stats count Organella chromosomes, sanity check here
    core_chr = core["coord_system"].get("chromosome", 0)
    core_adjusted_chrs = 0
    if core_chr:
        core_adjusted_chrs = core_chr - core_num_organella

    # Number of scaffolds from our core
    core_num_scaffolds = core["coord_system"].get("scaffold", 0)

    # NCBI includes the chromodomes in its stats
    core_adjusted_scaffolds = core_num_scaffolds + core_num_chrs

    # Compile the counts
    ncbi_counts = {
        "num_organella": len(ncbi_organella),
        "num_chromosomes": ncbi_main.get("total_number_of_chromosomes", 0),
        "num_scaffolds": ncbi_main.get("number_of_scaffolds", 0),
        "num_contigs": ncbi_main.get("number_of_contigs", 0),
    }
    core_counts = {
        "num_organella": core_num_organella,
        "num_chromosomes": core_adjusted_chrs,
        "num_scaffolds": core_adjusted_scaffolds,
        "num_contigs": core["coord_system"].get("contig", 0),
    }

    # Only compare contigs if there are any in the core db
    if core_counts["num_contigs"] == 0:
        del ncbi_counts["num_contigs"]
        del core_counts["num_contigs"]

    return _diff_dicts(ncbi_counts, core_counts)


def compare_annotation(ncbi: Dict, core: Dict) -> Dict:
    # Prepare counts to be comparable
    core_biotypes = core.get("genes", {}).get("biotypes", {})

    # We want a value for:
    # protein_coding = same
    # pseudogene = same
    # non_coding = total - rest
    # other = ?

    ncbi_counts = {
        "protein_coding": ncbi.get("protein_coding", 0),
        "pseudogene": ncbi.get("pseudogene", 0),
        "total_genes": ncbi.get("total", 0),
        "other": ncbi.get("other", 0),
    }

    # Add all pseudogenes
    num_pseudogenes = 0
    for name, num in core_biotypes.items():
        if re.match(".*pseudogen.*", name):
            num_pseudogenes += num

    core_counts = {
        "protein_coding": core_biotypes.get("protein_coding", 0),
        "pseudogene": num_pseudogenes,
        "total_genes": core.get("genes", {}).get("total", 0),
        "other": 0,
    }

    return _diff_dicts(ncbi_counts, core_counts)


def compare_stats(ncbi: Dict, core: Dict) -> Dict:
    """Compare stats from NCBI and our core."""

    ncbi_annotation_stats = ncbi.get("annotation_info", {}).get("stats", {}).get("gene_counts", {})
    core_assembly_stats = core.get("assembly_stats")
    core_annotation_stats = core.get("annotation_stats")

    comp = {
        "assembly_diff": compare_assembly(ncbi, core_assembly_stats),
    }
    if core_annotation_stats is not None and ncbi_annotation_stats is not None:
        comp["annotation_diff"] = compare_annotation(ncbi_annotation_stats, core_annotation_stats)

    return comp


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
        try:
            ncbi_stats = json.load(ncbi_fh)["reports"][0]
        except KeyError:
            ncbi_stats = {}
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
