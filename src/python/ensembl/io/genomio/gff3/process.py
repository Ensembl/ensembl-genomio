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
"""Simplify and fix a GFF3 file and returns both a cleaned up GFF3 file and a functional annotation
JSON file.
"""

import logging
from pathlib import Path

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

from .simplifier import GFFSimplifier
from .gene_merger import GFFGeneMerger


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description=(
            "Standardize the gene model representation of a GFF3 file, and extract the functional "
            "annotation in a separate file."
        )
    )
    parser.add_argument_src_path("--in_gff_path", required=True, help="Input GFF3 file")
    parser.add_argument_src_path("--genome_data", required=True, help="Genome JSON file")
    parser.add_argument(
        "--fail_missing_stable_ids", action="store_true", help="Do not generate IDs when missing/invalid"
    )
    parser.add_argument_dst_path("--out_gff_path", default=Path("gene_models.gff3"), help="Output GFF3 file")
    parser.add_argument_dst_path(
        "--out_func_path",
        default=Path("functional_annotation.json"),
        help="Output functional annotation JSON file",
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    # Merge multiline gene features in a separate file
    logging.info("Checking for genes to merge...")
    interim_gff_path = Path(f"{args.in_gff_path}_INTERIM_MERGE")
    merger = GFFGeneMerger()
    merged_genes = merger.merge(args.in_gff_path, interim_gff_path)
    num_merged_genes = len(merged_genes)
    in_gff_path = args.in_gff_path
    # If there are split genes, decide to merge, or just die
    if num_merged_genes > 0:
        # Report the list of merged genes in case something does not look right
        logging.info(f"{num_merged_genes} genes merged")
        logging.debug("\n".join(merged_genes))
        # Use the GFF with the merged genes for the next part
        in_gff_path = interim_gff_path

    # Load GFF3 data and write a simpler version that follows our specifications as well as a
    # functional annotation JSON file
    logging.info("Simplify and fix GFF3")
    gff_data = GFFSimplifier(args.genome_data)
    if args.fail_missing_stable_ids:
        gff_data.stable_ids.make_missing_stable_ids = False
    gff_data.simpler_gff3(in_gff_path)
    gff_data.records.to_gff(args.out_gff_path)
    gff_data.annotations.to_json(args.out_func_path)


if __name__ == "__main__":
    main()
