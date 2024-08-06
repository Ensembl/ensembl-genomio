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


def summary_feature_stats(in_file):
    """Run analysis of GFF file and produce a summary of feature types"""

    examiner = GFFExaminer()
    in_handle = open(in_file)

    pprint.pprint(examiner.available_limits(in_handle))
    print("\n\n")
    in_handle.close()
    sys.exit(0)

## New function to do main overlap detection
# def gff_digest()->None: #UPDATE RETURN TYPE


def scan_tree(intervals):
    """construct an interval tree using supplied genomic intervals, check all elements on the tree against iself and return any that hit 2 or more intervals (i.e. itself + 1 other)"""

    retlist = set()
    t = IntervalTree(Interval(*iv) for iv in intervals)

    for g in intervals:
        if len(t.overlap(g[0], g[1])) > 1:
            #            print( t.overlap( g[0], g[1]) )
            o = t.overlap(g[0], g[1])
            for x in o:
                retlist.add(x.data)

    return retlist


def write_report(out_file, seqDict, genesDict):
    """write final report"""
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

    except:
        logging.critical("failed to write to output file")

    return overlap_count


def get_intervals(record, genesDict, seqDict, seqname):
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
            seqDict[seqname]["plus"].append((int(r.location.start), int(r.location.end), str(r.id)))

        elif r.location.strand == -1:
            seqDict[seqname]["minus"].append((int(r.location.start), int(r.location.end), str(r.id)))

        else:
            msg = f"something went horribly wrong with the strand processing\n"
            logging.critical(msg)
            print(msg)

def gff_to_overlaps(input_GFF, parse_option, outfile, filter) -> None:
    """Master function to process input GFF file.
        
        Args:
            input_GFF:
            parse_option:
            outfile:
            filter:
        
        Returns:
            None
    """
    ## MAIN EDITS FOLLOW:
    logging.info("start")
    logging.info("gff file = " + input_GFF)
    logging.info("output file = " + outfile)
     
    if parse_stats_only is not None:
        logging.info("Alt processing: Not parsing the GFF, producing summary feature stats instead!")
        # Now do stats instead:
        summary_feature_stats(input_GFF)
    else:
        logging.info("Processing sequence feature overlaps!")
        gff_digest(args.input_GFF)


    
    ## remaining call on main edits:

    in_handle = open(in_file)
    limit_info = dict(gff_type=[filter])

    seqDict = defaultdict(dict)
    genesDict = {}

    for rec in GFF.parse(in_handle, limit_info=limit_info):
        #        print( type(rec) )
        #        print(rec.id)
        seqname = str(rec.id)
        if seqname not in seqDict:
            seqDict[seqname]["plus"] = []
            seqDict[seqname]["minus"] = []

        get_intervals(rec, genesDict, seqDict, seqname)

    in_handle.close()

    overlaps = write_report(out_file, seqDict, genesDict)

    msg_total_genes = f"total genes scanned = { len(genesDict) }"
    print(msg_total_genes)
    logging.info(msg_total_genes)

    msg_total_overlaps = f"total overlaps detected = { overlaps }"
    print(msg_total_overlaps)
    logging.info(msg_total_overlaps)

    logging.info("stop")
    sys.exit(0)




def main(in_file, out_file, filter):
    """Module entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_gff", 
        type=Path,
        required=True,
        default=None,
        help="Path of gff file to scan for overlaps")
    # processing_group = parser.add_mutually_exclusive_group(required=True)
    # processing_group.add_argument(
    #     "--gff_overlaps", 
    #     type=bool,
    #     required=False,
    #     action='store_true',
    #     help="Fully parse GFF and locate seqFeature overlaps.")
    # processing_group.add_argument(
    parser.add_argument(
        "--stats_only", 
        type=bool,
        required=False,
        action="store_true",
        help="Provide summary of GFF feature types but do not parse any GFF overlaps.")
    parser.add_argument(
        "--output",
        type=Path,
        required=False,
        default="seq_feature_overlaps.txt",
        help="path of file for output (default=overlaps)",
    )
    parser.add_argument(
        "--isolate_feature",
        type=str,
        required=False,
        default="gene",
        help="The seqfeature type to use for overlap isolation (default=gene)",
    )

    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()

    # Check optional processing param
    if args.stats_only is False:
        processing_selection = "detect_overlaps"
    else:
        logging.info("Summarizing file instead")
        processing_selection = "summarize_features"


    # Master GFF parsing function
    gff_to_overlaps(args.input_gff, processing_selection, args.output, args.isolate_feature)
