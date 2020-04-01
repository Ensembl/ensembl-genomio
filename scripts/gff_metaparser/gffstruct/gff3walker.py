import sys

from BCBio import GFF

class GFF3Walker:
  # BaseWalker?
  # single gff walking code
  # flag to iterate through qualifiers
  # tag=fullPath|leafQual|leafOrQual # class?? enum? named tuple?
  # /gene/mrna/exon
  # exon/id
  # exon
  def __init__(in_file, structure_tags="full", global_ctx = None):
    self.in_file = in_file
    self.gctx = global_ctx
    pass

  def walk(self, parser, out_file):
    gff =  GFF.parse(self.in_file)
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
            json_out.append(gene_info)
      GFF.write([out_rec], out_file)
      pass

  def process_feature(self, ft, out, prev_levels = []):
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

