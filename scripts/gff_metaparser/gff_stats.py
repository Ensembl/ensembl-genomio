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


import argparse
import io
import json
import sys


# locals
from gffstruct.statskeeper import StatsKeeper
from gffstruct.gff3walker import GFF3Walker
from gffstruct.utils import SeqLenDict
from gffstruct.validstruct import ValidStructures


def get_args():
    parser = argparse.ArgumentParser()
    # configs
    parser.add_argument(
        "--conf",
        metavar="structures.conf",
        required=False,
        type=argparse.FileType("rt", encoding="UTF-8"),
        help="valid structures config file",
    )
    # various options
    parser.add_argument("--stats_only", action="store_true", required=False, help="produce only stats output")
    parser.add_argument(
        "--no_id_stats", action="store_true", required=False, help="do not dump stemmed id stats"
    )
    parser.add_argument(
        "--fail_unknown", action="store_true", required=False, help="fail if unknown structure met"
    )
    parser.add_argument(
        "--processed_qualifiers",
        metavar="source,name,parent,dbxref,phase,product,protein_id,biotype",
        required=False,
        type=str,
        default="",
        help="comma separated list of qualifiers to use [keeping all by default]",
    )
    # ? separate sources, value to include #ALL_SOURCES with #
    parser.add_argument(
        "--no_longest_pfx", action="store_true", required=False, help="no lognest pfrefixes for IDs"
    )
    parser.add_argument(
        "--dump_used_options", action="store_true", required=False, help="dump used (not None) options"
    )
    parser.add_argument(
        "--rule_options",
        metavar="load_pseudogene_with_CDS,option2,...",
        required=False,
        action="append",
        help="options to control conditional rules",
    )
    parser.add_argument(
        "--no_contig_len_extenstion",
        action="store_true",
        required=False,
        help="do not extend contig length based on the feature boundaries (try to omit, if fails)",
    )
    # output
    parser.add_argument(
        "--stats_out",
        metavar="stats.out",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        default=sys.stdout,
        help="stats output [STDOUT]",
    )
    parser.add_argument(
        "--detailed_report",
        metavar="detailed.log",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        help="file to report every failed and modified model",
    )
    parser.add_argument(
        "--gff_out",
        metavar="out.gff3",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        default=sys.stdout,
        help="resulting gff output [STDOUT]",
    )
    # input
    parser.add_argument(
        "--fasta",
        metavar="fasta.fna",
        required=False,
        type=str,
        help="fasta file with sequences to calculate length (region length or max feature end will be used, if absent)",
    )
    parser.add_argument(
        "gff_in",
        metavar="in.gff3",
        type=argparse.FileType("rt", encoding="UTF-8"),
        help="input gff file (use '-' to read from STDIN)",
    )
    args = parser.parse_args()
    if args.dump_used_options:
        print("used options:", file=sys.stderr)
        for o, v in args.__dict__.items():
            if v:
                print("  %s: %s" % (o, v), file=sys.stderr)
    return args


def main():
    args = get_args()

    rule_options = args.rule_options
    if rule_options:
        rule_options = list(filter(None, map(lambda s: s.strip(), ",".join(rule_options).split(","))))

    parser = ValidStructures(args.conf, rule_options=rule_options)

    seq_len = SeqLenDict(args.fasta)

    no_id_stats_rules = ["UNSEEN", "IGNORE"]  # no id stats for these rules
    stats_keeper = StatsKeeper(
        detailed=args.detailed_report, no_longest_pfx=args.no_longest_pfx, no_id_stats_rules=no_id_stats_rules
    )
    gff3_walker = GFF3Walker(parser, args.gff_in, structure_tags="fullPath", global_ctx=stats_keeper)

    gff_out = args.gff_out
    if args.stats_only:
        gff_out = None
    gff3_walker.walk(
        out_file=gff_out, seq_len_dict=seq_len, contig_length_extension=not args.no_contig_len_extenstion
    )

    stats_keeper.dump(args.stats_out, id_stats=not args.no_id_stats)

    if args.fail_unknown:
        unseen_stats = stats_keeper.summary("UNSEEN")
        if unseen_stats:
            print("\nfailing because of the UNSEEN structure(s):\n%s" % unseen_stats, file=sys.stderr)
            sys.exit(1)


# main
if __name__ == "__main__":
    # execute only if run as a script
    main()

# zcat data/pfal/Pfalciparum.gff.gz | python new-genome-loader/scripts/gff_metaparser/gff_stats.py  --processed_qualifiers source,name,parent,dbxref,phase,product,protein_id,biotype --conf new-genome-loader/config/gff_metaparser/valid_structures.conf -

# zcat data/pfal/Pfalciparum.gff.gz | python new-genome-loader/scripts/gff_metaparser/gff_stats.py  --stats-only --processed_qualifiers source,name,parent,dbxref,phase,product,protein_id,biotype -
