import sys

class FixAction:
  def __init__(self, raw, rule):
     self._raw = raw
     self._rule_name = rule.NAME
     self._rule_lineno = rule._lineno
     self._action = self.parse(self._raw)

  def __bool__(self):
    return self._action is not None

  def __repr__(self):
    if self._action is None:
      return str(None)
    return "action: %s of %s at %s" % (
      self._raw, self._rule_name, self._rule_lineno
    )

  def parse(self, raw):
    out = []
    raw_splitted = raw.split("/")
    for p in raw.split("/"):
      if p.strip() == "" or p == "-":
        out.append({ "action": "exclude" })
        continue
      action = "rename"
      if p[0] == "+":
        p = p[1:]
        action = "add"
      _type, *quals = p.replace(",",".").split(".")
      _type = _type.strip()
      if _type == "":
        print("empty type for %s" % (self), file=sys.stderr)
        return None
      quals = dict([ (kv+"=").split("=")[0:2] for kv in quals if kv ])
      out.append({
        "action" : action,
        "type" : _type,
        "quals" : quals,
      })
    # TODO balanced number of all - added == number of old
    return out or None

  def act(self, ctx):
    re_ctx = ctx.get("_RECTX") and ctx["_RECTX"].groupdict() or None
    print("run_postponed for sub rectx:", re_ctx, ctx.get("_FULLTAG"), self._action, file=sys.stderr)
    return {}
    # run_postponed for sub rectx: {'MRNA': 'mRNA', 'CDS': 'CDS'} gene/mRNA/CDS ['gene/transcript.biotype=@MRNA/@CDS']

    #gene/@MRNA/@CDS	SUB	gene/transcript.biotype=@MRNA/@CDS
    #
    #gene/primary_transcript/mirna/exon	SUB	ncRNA_gene/pre_miRNA/-/exon
    #gene/primary_transcript/mirna/exon	SUB	ncRNA_gene/-/miRNA/exon
    #mirna	SUB	+ncRNA_gene/miRNA.biotype=miRNA/+exon
    #


