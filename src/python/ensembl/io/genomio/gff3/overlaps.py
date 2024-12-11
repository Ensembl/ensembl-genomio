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
"""Scan a GFF3 file to detect overlapping SeqFeature objects. Default object level => gene."""

__all__ = [
    "summarize_feature_stats",
    "identify_feature_overlaps",
    "scan_tree",
    "get_intervals",
]

from collections import defaultdict
import json
import logging
from pathlib import Path
from pprint import pprint

from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree

import ensembl.io.genomio
from ensembl.io.genomio.utils import print_json
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def summarize_feature_stats(gff_in: Path) -> None:
    """Analyse a GFF3 file and produce a summary of its feature types.

    Args:
        gff_in: User supplied GFF3 input file.
    """

    logging.info("Alt processing: Not parsing the GFF3, producing summary feature stats instead!")

    examiner = GFFExaminer()
    with gff_in.open("r", encoding="utf-8") as input_handle:
        pprint(examiner.available_limits(input_handle))
    input_handle.close()


def identify_feature_overlaps(gff_in: Path, output_file: Path, isolate_feature: str) -> None:
    """Detect overlapping GFF3 SeqFeature objects and dump to a report.

    Args:
        gff_in: User supplied GFF3 input file.
        output_file: Output file to write feature overlaps.
        isolate_feature: Sequence feature type to filter by.
    """
    logging.info("Processing sequence feature overlaps!")
    logging.info(f"Output file = {str(output_file)}")
    logging.info(f"Features filtered by type: {isolate_feature}")

    gff_type_filter: dict = {"gff_type": [isolate_feature]}
    seq_dict: dict = defaultdict(dict)
    genes_dict: dict = {}
    with gff_in.open("r", encoding="utf-8") as input_handle:
        for record in GFF.parse(input_handle, limit_info=gff_type_filter):
            seq_name = str(record.id)
            if seq_name not in seq_dict:
                seq_dict[seq_name]["plus"] = []
                seq_dict[seq_name]["minus"] = []

            get_intervals(record, genes_dict, seq_dict, seq_name)

    overlap_count = _write_report(output_file, seq_dict, genes_dict)

    result_total_features = f"In total {len(genes_dict)} {isolate_feature} features were scanned."
    print(result_total_features)
    logging.info(result_total_features)

    result_total_overlaps = f"In total {overlap_count} overlaps were detected."
    print(result_total_overlaps)
    logging.info(result_total_overlaps)

    logging.info("Finished all processing.")


def scan_tree(feature_intervals: list) -> set:
    """Construct an interval tree using supplied genomic intervals, check all elements on the tree against
    itself and return any that hit 2 or more intervals (i.e. itself + 1 other)

    Args:
        feature_intervals: Genome features to examine for coordinate (start/end) overlaps.

    Return:
        Set of intervals identified in the input GFF3 file that overlap with 2 or more intervals.
    """

    interval_sets = set()
    traversed_tree = IntervalTree(Interval(*iv) for iv in feature_intervals)

    for interval in feature_intervals:
        if len(traversed_tree.overlap(interval[0], interval[1])) > 1:
            overlap_interval = traversed_tree.overlap(interval[0], interval[1])

            for features in overlap_interval:
                interval_sets.add(features.data)

    return interval_sets


def _write_report(out_file: Path, seq_dict: dict, genes_dict: dict) -> int:
    """Write the final overlap report to output file.

    Args:
        out_file: Output name of file to dump detected feature overlaps.
        seq_dict: Sequence features.
        genes_dict: Unique (+ | - stranded) overlap features.

    Returns:
        Count of overlaps detected
    """
    overlap_count = 0
    overlap_features = []

    for seq_name in seq_dict:
        logging.info(f"{seq_name} plus  {len(seq_dict[seq_name]['plus'])}")
        logging.info(f"{seq_name} minus {len(seq_dict[seq_name]['minus'])}")

        positive_hit = scan_tree(seq_dict[seq_name]["plus"])
        logging.info(f"{len(positive_hit)} positive strand overlaps detected")

        negative_hit = scan_tree(seq_dict[seq_name]["minus"])
        logging.info(f"{len(negative_hit)} negative strand overlaps detected")

        uniq_features: set = positive_hit.union(negative_hit)
        overlap_count = overlap_count + len(uniq_features)

        for feature in uniq_features:
            format_feature = str(genes_dict[feature]).replace("'", '"')
            overlap_features.append(json.loads(format_feature))
        print_json(out_file, overlap_features)

    return overlap_count


def get_intervals(record: SeqRecord, genes_dict: dict, seq_dict: dict, seq_name: str) -> None:
    """Extract start/stop feature coordinates for use in creating intervaltree object.

    Args:
        record: Individual sequence record.
        genes_dict: Genes.
        seq_dict: Sequences.
        seq_name: Feature sequence name.
    """

    for feature in record.features:
        genes_dict[str(feature.id)] = {
            "sequence": f"{record.id}",
            "start": f"{int(feature.location.start) + 1}",
            "end": f"{int(feature.location.end)}",
            "strand": f"{feature.location.strand}",
            "name": f"{feature.id}",
        }

        if feature.location.strand == 1:
            seq_dict[seq_name]["plus"].append(
                (int(feature.location.start), int(feature.location.end), str(feature.id))
            )
        elif feature.location.strand == -1:
            seq_dict[seq_name]["minus"].append(
                (int(feature.location.start), int(feature.location.end), str(feature.id))
            )
        else:
            logging.critical("Something went wrong with the strand processing!")


def main() -> None:
    """Module entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    # Create parser with common arguments to be used by both subparsers
    base_parser = ArgumentParser(add_help=False)
    base_parser.add_argument_src_path("--input_gff", required=True, help="path of GFF3 file to process")
    base_parser.add_log_arguments(add_log_file=True)
    # Add subparsers with their parent being the base parser with the common arguments
    subparsers = parser.add_subparsers(title="Parse GFF3 and ", required=True, dest="subcommand")
    _ = subparsers.add_parser("stats", parents=[base_parser], help="Provide summary of feature types")
    overlaps_parser = subparsers.add_parser("overlaps", parents=[base_parser], help="Find feature overlaps")
    overlaps_parser.add_argument_dst_path(
        "--output_file", default="feature_overlaps.txt", help="path of output file"
    )
    overlaps_parser.add_argument(
        "--filter_type", default="gene", help="sequence feature type used for overlap isolation"
    )

    args = parser.parse_args()
    init_logging_with_args(args)

    logging.info("Starting processing...")
    logging.info(f"GFF input file = {str(args.input_gff)}")

    # Check optional processing param
    if args.subcommand == "stats":
        summarize_feature_stats(args.input_gff)
    else:
        identify_feature_overlaps(args.input_gff, args.output_file, args.filter_type)


if __name__ == "__main__":
    main()
