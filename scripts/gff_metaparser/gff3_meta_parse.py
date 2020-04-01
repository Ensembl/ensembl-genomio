import argparse
import io
import json
import sys

from BCBio import GFF
from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

# locals
from gffstruct.fannkeeper import FannKeeper
from gffstruct.gff3walker import ExtMapper
from gffstruct.gff3walker import IdNormalizer
from gffstruct.gff3walker import GFF3Walker
from gffstruct.metaparserstruct import MetaParserStructures


def get_args():
  parser = argparse.ArgumentParser()
  # various configs and maps
  parser.add_argument("--conf", metavar="gff_metaparser.conf", required=True,
                      type=argparse.FileType('rt', encoding='UTF-8'),
                      help="valid parser config file")
  parser.add_argument("--pfx_trims", metavar="gene:gene-,mrna:rna-gnl|WGS:VCGU|,exon:exon-gnl|WGS:VCGU|TCAL_,cds:cds",
                      required = False, type=str,
                      help="list of feature:id_pfx(,...) to trim id_pfx from feature's ID (reuse feature for multiple pfx, case ignored")
  parser.add_argument("--xref_map", metavar="xref_map.tsv", required=False,
                      type=argparse.FileType('rt', encoding='UTF-8'),
                      help="tab separated file with xref mappings")
  parser.add_argument("--xref_map_str", metavar="NCBI_GP:GenBank,gff_xref:ensembl_xref", required=False,
                      type = str,
                      help="comma separated list of xref mappings in the format from:to")
  # output
  parser.add_argument("--gff_out", metavar="out.gff3", required = False,
                      type=argparse.FileType('w', encoding='UTF-8'), default=sys.stdout,
                      help="resulting gff output [STDOUT]" )
  parser.add_argument("--fann_out", metavar="functional_annotation.json", required = False,
                      type=argparse.FileType('w',  encoding='UTF-8'), default=sys.stdout,
                      help="stats output [STDOUT]" )
  # input
  parser.add_argument("--fasta", metavar="fasta.fna", required = False,
                      type=str,
                      help="fasta file with sequences to calculate length (region length or max feature end will be used, if absent)")
  parser.add_argument("gff_in", metavar="in.gff3",
                      type=argparse.FileType('rt', encoding='UTF-8'),
                      help="input gff file (use '-' to read from STDIN)" )
  #
  args = parser.parse_args()
  return args




class PfxTrimmer(IdNormalizer):
  def __init__(self, trim_str):
    # features, can be reused, gen array
    pass


def main():
  args = get_args()

  parser = MetaParserStructures(args.conf)
  pfx_trimmer = PfxTrimmer(args.pfx_trims)
  xref_map = ExtMapper("xref", map_file = args.xref_map, map_str=args.xref_map_str)

  fann_ctx = FannKeeper()
  gff3_walker = GFF3Walker(args.gff_in,
    structure_tags = "leafQual", global_ctx = fann_ctx, norm_id = pfx_trimmer)

  gff3_walker.walk(parser)

  gff3_walker.dump_res_gff3(args.gff_out, maps = {"@xrefs": [xref_map]}) # or xref mapper obj
  fann_ctx.dump_json(args.fann_out)


  clean_tags = frozenset(["ID", "Parent"])
  def name_clean(x, tag = None):
    if tag and tag in clean_tags:
      return x.replace("cds-", "").replace("gene-", "").replace("rna-gnl|WGS:VCGU|", "").replace("exon-gnl|WGS:VCGU|", "")
    return x


if __name__ == "__main__":
    # execute only if run as a script
    main()


# cat tcal.gff3  | python new-genome-loader/scripts/gff_metaparser/gff3_tcal_parse.py --fann_out t.json -  > t.gff3
