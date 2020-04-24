import sys
from collections import defaultdict

from BCBio import GFF

from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


class SeqFeatComposer:
  # storing/building gff structures
  def __init__(self):
    self._data = defaultdict(lambda: defaultdict(dict))

  def gff_write(self, context, contig_len = None, out_file = None):
    if out_file is None:
      return
    out_rec_features = []
    self.compose(context, out_rec_features)
    if out_rec_features:
      out_rec = SeqRecord(UnknownSeq(length=self._contig_len), id = self._contig_id)
      out_rec.features = out_rec_features
      GFF.write([out_rec], out_file)

  def compose(self, context, out = None):
    if out is None:
      return
    for ctx in context.used_leaves():
      # self.merge_used_quals(ctx)
      print("leaf", list(map(lambda t: ctx.get(t), ["_ID", "_TYPE", "_FULLTAG", "_DEPTH"])), file = sys.stderr)

      used_quals = ctx.get("_GFF_USEDQUALS")
      print("leaf used_quals: ", used_quals, file = sys.stderr)

      parent_ctx = ctx.get("_PARENTCTX")
      while parent_ctx:
        print("parent", list(map(lambda t: parent_ctx.get(t), ["_ID", "_TYPE", "_FULLTAG", "_DEPTH"])), file = sys.stderr)
        used_quals = parent_ctx.get("_GFF_USEDQUALS")
        if (used_quals):
          print(used_quals, file = sys.stderr)
        parent_ctx = parent_ctx.get("_PARENTCTX")
      print("", file=sys.stderr)
    #
    # process fixes, posponed validations, etc, FIX ans SUB rules
    # replace out if fixes
    # update context / replace
    # gene/primary_transcript/mirna/exon
    # [gene , mrna utr exon exon cds cds  utr mirna exon , mrna utr exon cds utr, ...  ]
    # check that all the fixes are compatible, no different rules matched if  FIX and SUB rules (add rule method???)
    # store fixes into context
    # merges several fixes for different exons i.e.
    return
    if out:
      outft = SeqFeature(context["_LOCATION"],
        strand = context["_STRAND"],
        type = context["_TYPE"],
        qualifiers = context["_QUALS"])
      if res_sub_features:
        outft.sub_features = res_subfeatures
      out.append(outft)

