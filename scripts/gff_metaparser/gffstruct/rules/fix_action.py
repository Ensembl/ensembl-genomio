import sys

class FixAction:
  def __init__(self, raw, rule):
     self._raw = raw
     self._rule_name = rule.NAME
     self._rule_lineno = rule._lineno
     self._additions = 0
     self._exclusions = 0
     self._action = self.parse(self._raw)
     if self._additions > 0 and self._exclusions > 0:
       print("%s has add and exclude operations at the same time. skipping" % self, file=sys.stderr)
       self._action = None
     if not self._action:
       self._additions = 0
       self._exclusions = 0

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
        self._exclusions += 1
        continue
      action = "rename"
      if p[0] == "+":
        p = p[1:]
        action = "add"
        self._additions += 1
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
    return out or None

  def act(self, ctx):
    if self._action is None:
      return None

    fulltag = ctx.get("_FULLTAG")
    depth = ctx.get("_DEPTH")
    if depth != len(self._action) - self._additions:
      print("unbalanced number of tag %s  and actions %s. skipping", file = sys.stderr)
      return None

    re_ctx = ctx.get("_RECTX") and ctx["_RECTX"].groupdict() or None

    #print("run_postponed for sub rectx:", re_ctx, ctx.get("_FULLTAG"), self._action, file=sys.stderr)

    # run substitutions
    ait = self._action
    it = ctx
    for ait in reversed(self._action):
      aa = ait["action"]
      print(ait, file = sys.stderr)
      if aa == "add":
        continue
      if aa == "rename":
        _type = self.from_rectx(ait["type"], re_ctx)
        it["_TYPE"] = _type
        used_quals = it.get("_RULESDATA")["_ALL"].get("USEDQUALS")
        if used_quals:
          for q, v in ait.get("quals",{}).items():
            q = self.from_rectx(q, re_ctx)
            v = self.from_rectx(v, re_ctx)
            _q = q.lower()
            if _q in used_quals:
              if v is None or v == "":
                used_quals.pop(_q)
              else:
                tq, tv = used_quals[_q]
                used_quals[_q] = (tq, v)
            else:
              used_quals.update({_q:(q, v)})
          pass
      it = it.get("_PARENTCTX")
     # run delitions
     # mark as deleted ? add copy with skipped ?
    return {}

    #gene/@MRNA/@CDS	SUB	gene/transcript.biotype=@MRNA/@CDS
    #
    #gene/primary_transcript/mirna/exon	SUB	ncRNA_gene/pre_miRNA/-/exon
    #gene/primary_transcript/mirna/exon	SUB	ncRNA_gene/-/miRNA/exon
    #mirna	SUB	+ncRNA_gene/miRNA.biotype=miRNA/+exon
    #

  def from_rectx(self, x, re_ctx=None):
    if not re_ctx or x is None: return x
    if type(x) != str or len(x) < 2: return x
    if x[0] == "@":
      x = x[1:]
      x = re_ctx.get(x, x)
    return x


