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
"""Scan a GFF file to detect overlapping SeqFeature objects. Default object level => gene."""

__all__ = [
    "summarise_feature_stats",
    "gff_digest",
    "scan_tree",
    "write_report",
    "get_intervals",
    "gff_to_overlaps",
]

import logging
from pathlib import Path
import sys

from intervaltree import Interval, IntervalTree
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import pprint

import json  # ?? update
from collections import defaultdict

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def summarise_feature_stats(gff_in:Path) -> None:
    """Run analysis of GFF file and produce a summary of feature types
    
    Args:
        gff_in: 
    """

    examiner = GFFExaminer()
    with open(gff_in, "r") as input_handle:
        pprint.pprint(examiner.available_limits(input_handle))
    input_handle.close()


## New function to do main overlap detection
def gff_digest(gff_in, output_file:Path, filter:str) -> None:
    """
    Args:
        gff_in:
        output_file:

    Return:
        None
    """

    with open(gff_in, "r") as input_handle:
        limit_info = dict(gff_type=[filter])

        seqDict = defaultdict(dict)
        genesDict = {}

        for rec in GFF.parse(input_handle, limit_info=limit_info):
            seq_name = str(rec.id)
            if seq_name not in seqDict:
                seqDict[seq_name]["plus"] = []
                seqDict[seq_name]["minus"] = []

            get_intervals(rec, genesDict, seqDict, seq_name)

        input_handle.close()

        overlaps = write_report(output_file, seqDict, genesDict)

        msg_total_genes = (f"total genes scanned = { len(genesDict) }")
        print(msg_total_genes)
        logging.info(msg_total_genes)

        msg_total_overlaps = (f"total overlaps detected = { overlaps }")
        print(msg_total_overlaps)
        logging.info(msg_total_overlaps)

        logging.info("Finished all processing.")


def scan_tree(intervals)-> set:
    """Construct an interval tree using supplied genomic intervals, check all elements on the tree against 
    itself and return any that hit 2 or more intervals (i.e. itself + 1 other)
    
    Args:
        intervals:

    Returns:
        Set of intervals identified in the input GFF that overlaps with 2 or more intervals.
    
    """

    interval_sets = set()
    t = IntervalTree(Interval(*iv) for iv in intervals)

    for g in intervals:
        if len(t.overlap(g[0], g[1])) > 1:
            o = t.overlap(g[0], g[1])
            for x in o:
                interval_sets.add(x.data)

    return interval_sets


def write_report(out_file, seqDict, genesDict) -> None:
    """Write the final overlap report to output file.
    
    Args:
        out_file:
        seqDict:
        genesDict:

    Returns:
        None
    """
    overlap_count = 0

    try:
        with open(out_file, "w+") as fh:
            for x in seqDict:
                logging.info(f"{x} plus  { len( seqDict[x]['plus']) }")
                logging.info(f"{x} minus { len( seqDict[x]['minus']) }")

                pfound = scan_tree(seqDict[x]["plus"])
                logging.info(f"{ len(pfound) } positive strand overlaps detected")

                nfound = scan_tree(seqDict[x]["minus"])
                logging.info(f"{ len(nfound) } negative strand overlaps detected")

                uniq = pfound.union(nfound)
                overlap_count = overlap_count + len(uniq)

                for f in uniq:
                    fh.write(json.dumps(genesDict[f]) + "\n")
            fh.close()
    except IOError:
        logging.critical("failed to write to output file")

    return overlap_count


def get_intervals(record, genesDict, seqDict, seq_name) -> None:
    """extract start/stop feature coordinates for use in creating intervaltree object"""

    for r in record.features:
        genesDict[str(r.id)] = {
            "sequence": record.id,
            "start": int(r.location.start) + 1,
            "end": int(r.location.end),
            "strand": r.location.strand,
            "name": r.id,
        }

        if r.location.strand == 1:
            seqDict[seq_name]["plus"].append((int(r.location.start), int(r.location.end), str(r.id)))

        elif r.location.strand == -1:
            seqDict[seq_name]["minus"].append((int(r.location.start), int(r.location.end), str(r.id)))

        else:
            msg = "something went horribly wrong with the strand processing"
            logging.critical(msg)
            print(msg)


def gff_to_overlaps(input_GFF:Path, stats_flag:bool, output:Path, filter:str) -> None:
    """Master function to process input GFF file for overlaps or feature summary.

    Args:
        input_GFF:
        stats_flag:
        outfile:
        filter:

    Returns:
        None
    """
    
    ## MAIN EDITS FOLLOW:
    logging.info("Starting processing...")
    logging.info(f"GFF input file = {str(input_GFF)}")
    logging.info(f"Output file = {str(output)}")
    logging.info(f"Features filtered by type: {filter}")

    # Check optional processing param
    if stats_flag is True:
        logging.info("Alt processing: Not parsing the GFF, producing summary feature stats instead!")
        summarise_feature_stats(input_GFF)
    else:
        logging.info("Processing sequence feature overlaps!")
        gff_digest(input_GFF, output, filter)



def main() -> None:
    """Module entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_gff",
        type=Path,
        required=True,
        default=None,
        help="Path of gff file to scan for overlaps",
    )
    parser.add_argument(
        "--stats_only",
        required=False,
        action="store_true",
        help="Provide summary of GFF feature types but do not parse any GFF overlaps.",
    )
    parser.add_argument(
        "--output_file",
        type=Path,
        required=False,
        default="seqfeature_overlaps.txt",
        help="Path of file for output (default=overlaps).",
    )
    parser.add_argument(
        "--isolate_feature",
        type=str,
        required=False,
        default="gene",
        help="The seqfeature type used for overlap isolation.",
    )

    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    # Master GFF parsing function
    gff_to_overlaps(args.input_gff, args.stats_only, args.output_file, args.isolate_feature)
