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
"""Tool set to compare genome statistic between NCBI datasets and Ensembl's core databases."""

__all__ = ["stats_dict_cmp", "compare_assembly", "compare_annotation", "compare_stats", "compare_stats_files"]

import json
from os import PathLike
import re
from typing import Any, Dict

import ensembl.io.genomio
from ensembl.io.genomio.utils import get_json
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def stats_dict_cmp(ncbi: Dict[str, int], core: Dict[str, int]) -> Dict[str, Dict]:
    """Compares both dictionaries and returns the similar and different elements between both.

    The method assumes both dictionaries have the same set of keys. A key would be considered the
    same if its value in both dictionaries is the same, but will only be included in the returned
    dictionary if that value is different than 0.

    Args:
        ncbi: NCBI dataset statistics in key-value pairs.
        core: Core database statistics in key-value pairs.

    Returns:
        A dictionary with 2 keys:
        - "same": Pairs of key - value for those entries equal in both dictionaries.
        - "different": Keys that differ, with values for "ncbi", "core", and "diff", i.e. their
            difference represented as `core_value - ncbi_value`.

    """
    diff = {}
    same = {}
    for key, ncbi_count in ncbi.items():
        core_count = core[key]
        if ncbi_count == core_count:
            if ncbi_count != 0:
                same[key] = ncbi_count
        else:
            diff[key] = {"ncbi": ncbi_count, "core": core_count, "diff": core_count - ncbi_count}
    comparison: Dict[str, Dict] = {}
    if same:
        comparison["same"] = same
    if diff:
        comparison["different"] = diff
    return comparison


def compare_assembly(ncbi: Dict[str, Any], core: Dict[str, Any]) -> Dict[str, Dict]:
    """Extracts the assembly statistics and returns the comparison between both sources.

    The assembly statistics compared are the number of: organella, chromosomes, scaffolds and contigs.
    The last one is only included if NCBI's assembly is contig level.

    Args:
        ncbi: NCBI dataset assembly statistics.
        core: Core database assembly statistics.

    Returns:
        The common statistics with their value and the statistics with different value, including NCBI'
        and core database's values as well as their difference (`core_value - ncbi_value`).

    """
    # Prepare counts to be comparable to the NCBI stats
    ncbi_main = ncbi.get("assembly_stats", {})
    ncbi_info = ncbi.get("assembly_info", {})
    ncbi_organella = ncbi.get("organelle_info", [])

    # First count the organella
    core_num_organella = 0
    for loc, loc_count in core.get("locations", {}).items():
        if loc != "nuclear_chromosome":
            core_num_organella += loc_count

    # Our core stats count Organella chromosomes, sanity check here
    core_chr = core.get("coord_system", {}).get("chromosome", 0)
    core_adjusted_chrs = 0
    if core_chr:
        core_adjusted_chrs = core_chr - core_num_organella

    # Number of scaffolds from our core
    core_num_scaffolds = core.get("coord_system", {}).get("scaffold", 0)

    # NCBI includes the chromosomes in its scaffold count
    core_adjusted_scaffolds = core_num_scaffolds + core_adjusted_chrs

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
        "num_contigs": core.get("coord_system", {}).get("contig", 0),
    }

    # Only compare contigs if there are any in NCBI
    if ncbi_info.get("assembly_level") != "Contig":
        del ncbi_counts["num_contigs"]
        del core_counts["num_contigs"]

    return stats_dict_cmp(ncbi_counts, core_counts)


def compare_annotation(ncbi: Dict[str, Any], core: Dict[str, Any]) -> Dict[str, Dict]:
    """Extracts the annotation statistics and returns the comparison between both sources.

    Annotation statistics compared:
        - protein_coding
        - pseudogene (all pseudogene biotypes)
        - other (number of misc_RNA)
        - total

    Args:
        ncbi: NCBI dataset annotation statistics.
        core: Core database annotation statistics.

    Returns:
        The common statistics with their value and the statistics with different value, including NCBI'
        and core database's values as well as their difference (`core_value - ncbi_value`).

    """
    ncbi_counts = {
        "protein_coding": ncbi.get("protein_coding", 0),
        "pseudogene": ncbi.get("pseudogene", 0),
        "total_genes": ncbi.get("total", 0),
        "other": ncbi.get("other", 0),
    }

    # Prepare core database counts to be comparable
    core_biotypes = core.get("genes", {}).get("biotypes", {})

    # Add all pseudogenes
    num_pseudogenes = 0
    for name, num in core_biotypes.items():
        if re.match(".*pseudogen.*", name):
            num_pseudogenes += num

    # Other genes such as misc_mRNA
    num_others = core_biotypes.get("misc_RNA", 0)

    core_counts = {
        "protein_coding": core_biotypes.get("protein_coding", 0),
        "pseudogene": num_pseudogenes,
        "total_genes": core.get("genes", {}).get("total", 0),
        "other": num_others,
    }

    return stats_dict_cmp(ncbi_counts, core_counts)


def compare_stats(ncbi: Dict[str, Any], core: Dict[str, Any]) -> Dict[str, Dict]:
    """Compares the genome statistics between an NCBI dataset and a core database.

    Args:
        ncbi: NCBI dataset genome statistics.
        core: Core database genome statistics.

    Returns:
        The common statistics with their value and the statistics with different value, including NCBI'
        and core database's values as well as their difference (`core_value - ncbi_value`), for the
        assembly and annotation (if present in one of the sources) under "assembly_diff" and
        "annotation_diff" keys, respectively.

    """
    ncbi_annotation_stats = ncbi.get("annotation_info", {}).get("stats", {}).get("gene_counts", {})
    core_assembly_stats = core.get("assembly_stats", {})
    core_annotation_stats = core.get("annotation_stats", {})

    comp: Dict[str, Dict] = {
        "assembly_diff": compare_assembly(ncbi, core_assembly_stats),
    }
    if core_annotation_stats or ncbi_annotation_stats:
        comp["annotation_diff"] = compare_annotation(ncbi_annotation_stats, core_annotation_stats)
    return comp


def compare_stats_files(ncbi_file: PathLike, core_file: PathLike) -> Dict[str, Dict]:
    """Compares the genome statistics between an NCBI dataset and a core database.

    Args:
        ncbi_file: NCBI dataset genome statistics JSON file.
        core_file: Core database genome statistics JSON file.

    Returns:
        The common statistics with their value and the statistics with different value, including NCBI'
        and core database's values as well as their difference (`core_value - ncbi_value`), for the
        assembly and annotation (if present in one of the sources) under "assembly_diff" and
        "annotation_diff" keys, respectively.

    """
    ncbi_stats = {}
    ncbi_stats = get_json(ncbi_file)["reports"][0]
    core_stats = get_json(core_file)
    all_stats = compare_stats(ncbi_stats, core_stats)
    return all_stats


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Compares the genome statistics between an NCBI dataset and a core database."
    )
    parser.add_argument_src_path("--ncbi_stats", required=True, help="NCBI dataset stats JSON file")
    parser.add_argument_src_path("--core_stats", required=True, help="core database stats JSON file")
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()

    # Configure and initialise logging
    init_logging_with_args(args)

    report = compare_stats_files(args.ncbi_stats, args.core_stats)
    print(json.dumps(report, indent=2, sort_keys=True))
