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


from intervaltree import Interval, IntervalTree
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import pprint
import argparse
import sys
import logging
import json
from collections import defaultdict


def stats(in_file):
    """run analysis of GFF file and produce a summary of feature types"""

    examiner = GFFExaminer()
    in_handle = open(in_file)

    print(f"\nrunning analysis of GFF file\n")

    pprint.pprint(examiner.available_limits(in_handle))
    print("\n\n")
    in_handle.close()
    sys.exit(0)


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


def main(in_file, out_file, filter):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Scan a GFF file to detect overlapping seqfeature objects (default = gene)."
    )
    parser.add_argument("--gff", help="path of gff file to scan for overlaps", required=True, nargs="?")
    parser.add_argument(
        "--stats", help="run stats analysis of GFF file instead of parsing GFF file", action="store_true"
    )
    parser.add_argument(
        "--logfile", help="path of file for logging (default=log.overlaps)", nargs="?", default="log.overlaps"
    )
    parser.add_argument(
        "--output",
        help="path of file for output (default=overlaps)",
        required=True,
        nargs="?",
        default="overlaps",
    )
    parser.add_argument(
        "--filter",
        help="seqfeature type to use for overlap isolation (default=gene)",
        nargs="?",
        default="gene",
    )
    args = parser.parse_args()

    if args.stats:
        stats(args.gff)
        sys.exit(0)

    logging.basicConfig(filename=args.logfile, level=logging.DEBUG, format="%(asctime)s %(message)s")
    logging.info("start")
    logging.info("gff file = " + args.gff)
    if args.output:
        logging.info("output file = " + args.output)
    else:
        logging.info("no output file supplied - assuming STDOUT")

    main(args.gff, args.output, args.filter)
