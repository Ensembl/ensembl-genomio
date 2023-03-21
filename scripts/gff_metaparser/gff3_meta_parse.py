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
from gffstruct.fannkeeper import FannKeeper
from gffstruct.gff3walker import GFF3Walker
from gffstruct.metaparserstruct import MetaParserStructures
from gffstruct.utils import ExtMapper
from gffstruct.utils import PfxTrimmer
from gffstruct.utils import SeqLenDict


def get_args():
    parser = argparse.ArgumentParser()
    # various configs and maps
    parser.add_argument(
        "--conf",
        metavar="gff_metaparser.conf",
        required=True,
        type=argparse.FileType("rt", encoding="UTF-8"),
        help="valid parser config file",
    )
    parser.add_argument(
        "--conf_patch", metavar="conf.patch", required=False, type=str, help="config patch file"
    )
    parser.add_argument(
        "--pfx_trims",
        metavar="ANY!:.+\|,ANY:id-,ANY:gene-,ANY:rna-,ANY:mrna-,cds:cds-,exon:exon-",
        required=False,
        type=str,
        help="""Comma separated list of `feature:id_pfx` pairs to trim. `feature` case is ignored.
                              If `feature` part is ommited or 'ANY' is used, `id_pfx` will be treamed from any `feature`.
                              If '!' is NOT used before ':', `id_pfx` is escaped.
                              Same `feature` (also true for 'ANY' or empty) can be specified several times.""",
    )
    parser.add_argument(
        "--xref_map",
        metavar="xref_map.tsv",
        required=False,
        type=argparse.FileType("rt", encoding="UTF-8"),
        help="tab separated file with xref mappings",
    )
    parser.add_argument(
        "--xref_map_str",
        metavar="NCBI_GP:GenBank,gff_xref:ensembl_xref",
        required=False,
        type=str,
        help="comma separated list of xref mappings in the format from:to",
    )
    parser.add_argument(
        "--dump_used_options", action="store_true", required=False, help="dump used (not None) options"
    )
    parser.add_argument(
        "--no_contig_len_extenstion",
        action="store_true",
        required=False,
        help="do not extend contig length based on the feature boundaries (try to omit, if fails)",
    )
    # output
    parser.add_argument(
        "--gff_out",
        metavar="out.gff3",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        default=sys.stdout,
        help="resulting gff output [STDOUT]",
    )
    parser.add_argument(
        "--fann_out",
        metavar="functional_annotation.json",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        default=sys.stdout,
        help="functional annnotation output [STDOUT]",
    )
    parser.add_argument(
        "--seq_region_out",
        metavar="seq_region_raw.json",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        default=sys.stdout,
        help="seq_region metadata json output [STDOUT]",
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
    #
    args = parser.parse_args()
    if args.dump_used_options:
        print("used options:", file=sys.stderr)
        for o, v in args.__dict__.items():
            if v:
                print("  %s: %s" % (o, v), file=sys.stderr)
    return args


def seq_region_filter(x, prop="object_type", tag="seq_region"):
    # todo: get tag/"seq_region" from config
    return x and prop in x and x[prop] == tag or False


def main():
    args = get_args()

    parser = MetaParserStructures(args.conf, conf_patch=args.conf_patch)

    pfx_trimmer = PfxTrimmer(args.pfx_trims)
    seq_len = SeqLenDict(args.fasta)
    xref_map = ExtMapper("xref/dbname", map_file=args.xref_map, map_str=args.xref_map_str)

    fann_ctx = FannKeeper()
    gff3_walker = GFF3Walker(
        parser, args.gff_in, structure_tags="anyQual", global_ctx=fann_ctx, norm_id=pfx_trimmer
    )

    gff3_walker.walk(
        out_file=args.gff_out, seq_len_dict=seq_len, contig_length_extension=not args.no_contig_len_extenstion
    )

    fann_ctx.dump(args.fann_out, maps=[xref_map], dump_filter=lambda x: not seq_region_filter(x))
    fann_ctx.dump(args.seq_region_out, dump_filter=lambda x: seq_region_filter(x))


# main
if __name__ == "__main__":
    # execute only if beeing run as a script
    main()


# cat tcal.gff3 | python ./new-genome-loader/scripts/gff_metaparser/gff3_meta_parse.py --conf ./new-genome-loader/config/gff_metaparser/metaparser.conf --pfx_trims 'ANY!:.+\|,ANY:gene-,cds:cds-,exon:exon-' --conf_patch ./new-genome-loader/config/gff_metaparser/metaparser/xref2gene.patch --fann_out t.json --seq_region_out sr.json - > t.gff3
