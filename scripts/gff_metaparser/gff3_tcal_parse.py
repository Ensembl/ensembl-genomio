import argparse
import io
import json
import pprint
import sys

from BCBio import GFF
from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

# xamining  data/GCA_007210705.1_Tcal_SD_v2.1/GCA_007210705.1_Tcal_SD_v2.1_genomic.gff.gz
# gene/mRNA/exon  74199   Genbank {"gene": ["gene-TCAL_", 74199], "mRNA": ["rna-gnl|WGS:VCGU|TCAL_", 74199], "exon": ["exon-gnl|WGS:VCGU|TCAL_", 74199]}
# gene/mRNA/CDS   68827   Genbank {"gene": ["gene-TCAL_", 68827], "mRNA": ["rna-gnl|WGS:VCGU|TCAL_", 68827], "CDS": ["cds-TRY", 68827]}
#region  459     Genbank {"region": ["", 459]}

def get_args():
  parser = argparse.ArgumentParser()
  # ? separate sources, value to include #ALL_SOURCES with #
  # file arguments
  parser.add_argument("--fann_out", metavar="functional_annotation.json", required = False,
                      type=argparse.FileType('w',  encoding='UTF-8'), default=sys.stdout,
                      help="stats output [STDOUT]" )
  parser.add_argument("--gff_out", metavar="out.gff3", required = False,
                      type=argparse.FileType('w', encoding='UTF-8'), default=sys.stdout,
                      help="resulting gff output [STDOUT]" )
  parser.add_argument("gff_in", metavar="in.gff3",
                      type=argparse.FileType('rt', encoding='UTF-8'),
                      help="input gff file (use '-' to read from STDIN)" )
  args = parser.parse_args()
  return args



def main():
  args = get_args()

  interest_q = {
    "ID" : "ID",
    #"Parent" : "Parent",
    "gene_biotype" : "biotype",
    "phase" : "phase",
  }

  clean_tags = frozenset(["ID", "Parent"])
  def name_clean(x, tag = None):
    if tag and tag in clean_tags:
      return x.replace("cds-", "").replace("gene-", "").replace("rna-gnl|WGS:VCGU|", "").replace("exon-gnl|WGS:VCGU|", "")
    return x

  xrefs_map = {
    "NCBI_GP" : "GenBank",
  }

  json_type_map = {
    "gene" : "gene",
    "mrna" : "transcript",
  }
  # how to pass properties to upper object
  #  ? update prev_levels
  #  ? prev levels is for json ???


# {"object_type":"gene","id":"AGAP004678","xrefs":[{"info_type":"DIRECT","id":"agCG45746","dbname":"Celera_Gene"},{"info_type":"DIRECT","id":"ebi640","dbname":"Celera_Gene"},{"info_type":"DEPENDENT","id":"1275380","dbname":"EntrezGene"}]}


  def process_feature(ft, out, prev_levels = []):
    #quals_to_keep =
    #quals_to_json = 
    quals = {
      interest_q[k]: name_clean(v[0], k)
        for k, v in ft.qualifiers.items() if k in interest_q
    }

    if ft.type == "gene":
      prev_levels.append({"object_type": "gene", "id" : quals["ID"], "xrefs": []})

    if "Dbxref" in ft.qualifiers:
      # could be many
      xref = ft.qualifiers["Dbxref"][0]
      xsrc, xid = xref.split(":", 1)
      xsrc = xsrc in xrefs_map and xrefs_map[xsrc] or xsrc
      if prev_levels and xsrc:
        prev_levels[0]["xrefs"].append({
          "info_type":"DIRECT","id":xid,"dbname":xsrc,
        })

    outft = SeqFeature(ft.location, strand = ft.strand, type = ft.type, qualifiers = quals)

    outft.sub_features = []
    for sft in ft.sub_features:
      process_feature(sft, outft.sub_features, prev_levels)

    out.append(outft)

    return


  json_out = []
  gff =  GFF.parse(args.gff_in)
  for contig in gff:
    out_rec = SeqRecord(UnknownSeq(length=len(contig)), id = contig.id)
    out_rec.features = []
    for cnt, topft in enumerate(contig.features):
      if topft.type == "gene":
        info = []
        process_feature(topft, out_rec.features, info)
        if info:
          gene_info = info[0]
          if "xrefs" in gene_info:
            #gene_info["xrefs"] = list(set(gene_info["xrefs"]))
            if gene_info["xrefs"]:
              gene_info["xrefs"] = gene_info["xrefs"][0]
            if not gene_info["xrefs"]:
              del(gene_info["xrefs"])
          json_out.append([gene_info])
    GFF.write([out_rec], args.gff_out)

  if (json_out):
    json.dump(json_out, args.fann_out, indent = 2)


if __name__ == "__main__":
    # execute only if run as a script
    main()


# cat tcal.gff3  | python new-genome-loader/scripts/gff_metaparser/gff3_tcal_parse.py --fann_out t.json -  > t.gff3
