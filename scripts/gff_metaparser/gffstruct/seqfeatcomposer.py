import sys

from BCBio import GFF

from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


class SeqFeatComposer:
  # storing/building gff structures
  def __init__(self, gene_id_prepend):
    self._topf = set()
    self._processed = {}
    self._gene_id_prepend = gene_id_prepend

  def gff_write(self, out_rec_features, contig_id = None, contig_len = None, out_file = None):
    if out_file is None:
      return
    if out_rec_features:
      out_rec = SeqRecord(UnknownSeq(length=contig_len), id = contig_id)
      out_rec.features = out_rec_features
      GFF.write([out_rec], out_file)

  def cid(self, ctx):
    p_type, p_id, p_ctx = "", "", ctx.get("_PARENTCTX")
    if p_ctx:
      p_type, p_id = map(lambda k: p_ctx.get(k, ""), ["_TYPE", "_ID"])
    _type, _start, _end, _strand, _id = map(lambda k: ctx.get(k, ""), ["_TYPE", "_START", "_END", "_STRAND", "_ID"])
    return (p_type, p_id, _type, _start, _end, _strand, _id)

  def processed_add(self, ctx, used_quals = None, kid = None):
    cid = self.cid(ctx)
    feat = self._processed.get(cid)

    # flattern and exclude "parent"s
    used_quals_flat = { v[0]:v[1] for k, v in (used_quals or {}).items() if k != "parent" and k[0] != "_" }

    if not feat:
      location = ctx["_LOCATION"] # die if no location
      _type = ctx["_TYPE"] # should be already changed by rules
      is_leaf = ctx["_ISLEAF"]
      strand = ctx.get("_STRAND")
      phase = ctx.get("_PHASE")
      # gene:ID workaround
      update_id = False
      _id = ctx.get("_ID")
      if (_type == "gene" and self._gene_id_prepend):
        _id = "gene:%s" %_id
        update_id = True
      # quals
      quals = used_quals_flat or {}
      if not is_leaf and "ID" not in quals:
        update_id = True
      if update_id:
        quals["ID"] = _id
      if phase is not None and "phase" not in quals:
        quals["phase"] = phase
      # create feat object
      quals = self.sort_quals(quals)
      obj = SeqFeature(location, strand = strand, type = _type, qualifiers = quals) # NB referebce to quals / allows post assign modifications
      # fill processed
      feat = { "cid" : cid, "obj" : obj, "quals" : quals, "kids" : dict(), "is_leaf" : is_leaf }
      self._processed[cid] = feat
    else:
      # gene:ID workaround
      if "ID" in used_quals_flat:
        _id = used_quals_flat["ID"][0]
        if feat.get("obj") and feat["obj"].type == "gene" and self._gene_id_prepend:
          if not _id.startswith("gene:"):
            used_quals_flat["ID"][0] = "gene:%s" % _id
      # fill quals
      used_quals_flat = self.sort_quals(used_quals_flat)
      # assert( id(feat["quals"]) = id(feat["obj"].qualifiers) )
      # print(cid, hex(id(feat["quals"])), hex(id(feat["obj"].qualifiers)), file=sys.stderr)
      feat["quals"].update(used_quals_flat)

    if kid:
      kcid = kid["cid"]
      if kcid not in feat["kids"]:
        feat["kids"][kcid] = len(feat["kids"])

    return feat

  def sort_quals(self, quals):
    if not quals:
      return quals
    filtered_id = filter(lambda x: x[0] == "ID", quals.items())
    filtered = filter(lambda x: x[0][0] != "_" and x[0] != "ID", quals.items())
    return dict(list(filtered_id) + list(sorted(filtered, key=lambda x:x[0])))

  def compose(self, context, out = None):
    if out is None: return

    # clear
    self._topf = set()
    self._processed = {}

    # go from used leaves to top
    for ctx in context.used_leaves():
      # print("compose ctx type", ctx.get("_TYPE"), file =sys.stderr)
      rules_data = ctx.get("_RULESDATA")
      if rules_data and "_ALL" in rules_data and "USEDQUALS" in rules_data["_ALL"]:
        used_quals = rules_data["_ALL"]["USEDQUALS"]

      feat = self.processed_add(ctx, used_quals = used_quals)
      parent_ctx = ctx.get("_PARENTCTX")

      while parent_ctx:
        rules_data = parent_ctx.get("_RULESDATA")
        if rules_data and "_ALL" in rules_data and "USEDQUALS" in rules_data["_ALL"]:
          used_quals = rules_data["_ALL"]["USEDQUALS"]

        feat = self.processed_add(parent_ctx, used_quals = used_quals, kid = feat)
        parent_ctx = parent_ctx.get("_PARENTCTX")

      self._topf.add(feat["cid"])

    # update subfeatures
    for cid, feat in self._processed.items():
      if feat["kids"]:
        feat["obj"].sub_features = [ self._processed[cid]["obj"] for cid, _ in sorted(feat["kids"].items(), key = lambda k: -k[1]) ]
    # fill out
    for cid in sorted(self._topf):
      #print(cid, file = sys.stderr)
      #print(self._processed[cid], file = sys.stderr)
      #pass
      out.append(self._processed[cid]["obj"])

