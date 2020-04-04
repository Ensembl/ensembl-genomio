import re
import sys

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


class GFF3Walker:
  _valid_structure_tags = frozenset("fullPath leafQual".split()) # /gene/mrna/exon vs exon/id, exon

  def __init__(in_file, structure_tags="fullPath", global_ctx = None, norm_id = None):
    if structure_tags not in _valid_structure_tags:
      raise Exception("not a valid 'structure_tags' value: %s" % structure_tags)
    self._struct_tags = structure_tags
    self._in_file = in_file
    self._gctx = global_ctx
    self._norm_id = norm_id

  def norm_id(self, id, type):
    if self._norm_id is None:
      return id
    return self._norm_id(id, type = type)

  # parser = MetaParserStructures(args.conf)
  # gff3_walker = GFF3Walker(args.gff_in, structure_tags = "leafQual",
  #                           global_ctx = fann_ctx, norm_id = pfx_trimmer)
  # gff3_walker.walk(parser, out_file = args.gff_out, seq_len_dict = seq_len)
  def walk(self, parser, out_file = None, seq_len_dict = None):
    gff =  GFF.parse(self.in_file)
    if not out_file:
      pass

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


  def down_features(self, feature, prefix = None, prev_ids = []):
    #print(feature)
    # use https://docs.python.org/3/library/operator.html attrgetter methodcaller ??
    #print(feature.type, feature.id, feature.location)
    quals = feature.qualifiers
    out = { k : (len(v) > 0 and v[0] or None)
              for k,v in quals.items()
                if (not useful_qls or k.lower() in useful_qls) }
    #print(out)

    if prefix is None:
      prefix = feature.type
    else:
      prefix += "/" + feature.type
    new_prev_ids = prev_ids + [feature.id]

    if not feature.sub_features:
      known_structures.process(ValidStructures.Structure(prefix))

      # analyze
      stats[prefix] += 1
      stats_sources[prefix].add(out["source"])
      for (ftype, fid) in zip (prefix.split("/"), new_prev_ids):
        stats_prefixes[prefix][ftype].add(fid)

    for sf in feature.sub_features:
      down_features(sf, prefix, new_prev_ids)


  def get_stats(self):
    useful_qls = frozenset(filter(lambda x: x != "", args.only_qualifiers.split(",")))

    stats = defaultdict(int)
    stats_sources = defaultdict(set)
    stats_prefixes = defaultdict(lambda: defaultdict(PfxTr)) # should be paths:item based

    gff = GFF.parse(args.gff_in)
    for contig in gff:
      for cnt, feature in enumerate(contig.features):
        self.down_features(feature)
        # if cnt > 20: break
      # break


