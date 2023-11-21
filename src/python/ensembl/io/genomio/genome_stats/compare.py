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
"""Compare stats in a JSON from NCBI dataset and a JSON from our core db.

Returns the JSON from our core db including a section with comparisons.
"""

__all__ = ["compare_assembly", "compare_annotation", "compare_stats"]

import json
import re
from typing import Any, Dict
import logging

from ensembl.io.genomio.utils import get_json
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging


def _diff_dicts(ncbi: Dict[str, int], core: Dict[str, int]) -> Dict[str, Any]:
    """Compare two dicts with the same keys and compute the difference of their values.

    Returns:
        A dict with 2 dicts:
        - same: List of keys that have the same count
        - different: Dict of keys with different count, with values for ncbi, core, and diff (=core - ncbi)
    """
    diff = {}
    same = {}
    for key, ncbi_count in ncbi.items():
        core_count = core[key]
        if ncbi_count == core_count:
            if ncbi_count != 0:
                same[key] = ncbi_count
            continue
        diff[key] = {"ncbi": ncbi_count, "core": core_count, "diff": core_count - ncbi_count}

    comp: Dict[str, Any] = {}
    if same:
        comp["same"] = same
    if diff:
        comp["different"] = diff
    return comp


def compare_assembly(ncbi: Dict[str, Any], core: Dict[str, Any]) -> Dict[str, Any]:
    """Returns a compilation of count comparisons.
    Each comparison is a dict with the value from NCBI, from Core, and their diff.

    Args:
        ncbi: Dict of stats from NCBI.
        core: Dict of assembly stats from the core.

    Returns:
        Dict: Each count from NCBI, from Core, and their diff.
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

    # NCBI includes the chromosomes in its stats
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
    if ncbi_info.get("assembly_level") != "Contig":
        del ncbi_counts["num_contigs"]
        del core_counts["num_contigs"]

    return _diff_dicts(ncbi_counts, core_counts)


def compare_annotation(ncbi: Dict, core: Dict) -> Dict:
    """Compare NCBI vs Core annotation stats (biotype counts).

    Args:
        ncbi: Dict of biotype counts from NCBI.
        core: Dict of biotype counts from the core.

    Returns:
        Dict: Each count from NCBI, from Core, and their diff.
    """
    # Prepare counts to be comparable
    core_biotypes = core.get("genes", {}).get("biotypes", {})

    # We want a value for:
    # protein_coding = same
    # pseudogene = all pseudogene biotypes
    # total = same
    # other = number of misc_RNA

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

    # Others? misc_mRNA (or anything similar?)
    num_others = core_biotypes.get("misc_RNA", 0)

    core_counts = {
        "protein_coding": core_biotypes.get("protein_coding", 0),
        "pseudogene": num_pseudogenes,
        "total_genes": core.get("genes", {}).get("total", 0),
        "other": num_others,
    }

    return _diff_dicts(ncbi_counts, core_counts)


def compare_stats(ncbi: Dict, core: Dict) -> Dict[str, Any]:
    """Compare stats from NCBI and our core.

    Args:
        ncbi: Dict of stats from NCBI.
        core: Dict of stats from the core.

    Returns:
        Dict: Each count from NCBI, from Core, and their diff.
    """

    ncbi_annotation_stats = ncbi.get("annotation_info", {}).get("stats", {}).get("gene_counts", {})
    core_assembly_stats = core.get("assembly_stats", {})
    core_annotation_stats = core.get("annotation_stats", {})

    comp: Dict[str, Any] = {
        "assembly_diff": compare_assembly(ncbi, core_assembly_stats),
    }
    if core_annotation_stats and ncbi_annotation_stats:
        comp["annotation_diff"] = compare_annotation(ncbi_annotation_stats, core_annotation_stats)

    return comp


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Compare genome statistics between an NCBI dataset and a core database."
    )
    parser.add_argument_src_path("--ncbi_stats", required=True, help="NCBI dataset JSON file")
    parser.add_argument_src_path("--core_stats", required=True, help="Core database JSON file")
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()

    # Configure and initialise logging
    init_logging(args.log_level, args.log_file, args.log_file_level)

    try:
        ncbi_stats = get_json(args.ncbi_stats)["reports"][0]
    except (json.decoder.JSONDecodeError, KeyError):
        logging.warning(f"{args.ncbi_stats} JSON file is empty")
        ncbi_stats = {}
    else:
        logging.info(f"{args.ncbi_stats} JSON .'reports' obtained")
    core_stats = get_json(args.core_stats)
    all_stats = compare_stats(ncbi_stats, core_stats)

    print(json.dumps(all_stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
