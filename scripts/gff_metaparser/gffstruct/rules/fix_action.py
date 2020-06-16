import copy
import sys

from .valid import ValidRule

class FixAction:
  def __init__(self, raw, rule):
     self._raw = raw
     self._rule_name = rule.NAME
     self._rule_lineno = rule._lineno
     #
     self._additions = 0
     self._exclusions = 0
     self._copy_leaves = 0
     #
     self._action = None
     self._action = self.parse(self._raw)
     if self._additions > 0 and self._exclusions > 0:
       print("%s has add and exclude operations at the same time. skipping" % self, file=sys.stderr)
       self._action = None
     if not self._action:
       self._additions = 0
       self._exclusions = 0
       self._copy_leaves = 0

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
    for p in raw_splitted:
      if p.strip() == "" or p == "-":
        out.append({ "action": "exclude" })
        self._exclusions += 1
        continue
      action = "rename"
      if p[0] == "+":
        action = "add"
        p = p[1:]
        self._additions += 1
      elif p[0] == "!":
        action = "copy_leaf"
        if id(p) != id(raw_splitted[-1]):
          self._action = False
          print("copy_leaf(!) used not with leaf in %s" % (self), file=sys.stderr)
          return None
        p = p[1:]
        self._copy_leaves += 1
      _type, *quals = p.replace(",",".").split(".")
      _type = _type.strip()
      if _type == "":
        print("empty type for %s" % (self), file=sys.stderr)
        return None
      quals = dict([ (kv+"=").split("=")[0:2] for kv in quals if kv ])
      # remove duplicated id for the copied leaf
      if action == "copy_leaf": quals["id"] = None
      out.append({
        "action" : action,
        "type" : _type,
        "quals" : quals,
      })
    return out or None

  def act(self, ctx):
    """
      copy node if updating parentctx:
        if adding or changing parentctx (skipping)
        remove/set _ISLEAF when needed
      return dict added nodes {id(node):id}
      mark removed/updated leaves as {id(node):none}
      mark leaves as to be removed from old leaves
    """
    if self._action is None:
      return None
    fulltag = ctx.get("_FULLTAG")
    depth = ctx.get("_DEPTH")
    if depth != len(self._action) - self._additions:
      print("unbalanced number of tag %s  and actions %s. skipping", file = sys.stderr)
      return None

    re_ctx = ctx.get("_RECTX") and ctx["_RECTX"].groupdict() or None
    self.run_subsitutions(ctx, re_ctx)

    new_nodes = {}
    # add / remove, updating parent ctx
    prev, it = None, ctx
    for ait in reversed(self._action):
      is_leaf = it.get("_ISLEAF", False)
      parent = it.get("_PARENTCTX")
      aa = ait["action"]
      if aa == "copy_leaf":
        if is_leaf:
          it = self.copy_node(it, new_nodes, keep_leaf = True, clean = True)
          self.update_node(it, ait, re_ctx)
      elif aa == "exclude":
        if is_leaf:
          self.del_node_if_leaf(it, new_nodes)
          parent = self.copy_node(parent, new_nodes)
          parent["_ISLEAF"] = True
        elif prev:
          prev = self.copy_node(prev, new_nodes)
          prev["_PARENTCTX"] = parent
        #
        prev, it = prev, parent
      elif aa == "add":
        # use id sfx based on depth, mb different for gene
        pass
      prev, it = it, parent
    #
    return new_nodes

  def copy_node(self, node, new_nodes, keep_leaf = False, clean = False):
    if node.get("_ISCOPY"):
      return node
    ncopy = copy.deepcopy(node)
    ncopy["_ISCOPY"] = True
    new_nodes[id(ncopy)] = ncopy
    # clean
    if clean:
      pass
      ncopy.get("_RULESDATA")["_ALL"]["USEDQUALS"] = {}
    # mark old leaf as unused
    if not keep_leaf:
      self.del_node_if_leaf(node, new_nodes)
    return ncopy

  def del_node_if_leaf(self, node, new_nodes):
    if node.get("_ISLEAF", False):
      new_nodes[id(node)] = None

    #gene/@MRNA/@CDS	SUB	gene/-/@CDS
    # mirna SUB +ncRNA_gene/miRNA.biotype=miRNA/+exon
    # gene/pre_mirna SUB ncRNA_gene/pre_miRNA/+exon
    # @RNA_EXON_ONLY/exon SUB +gene.biotype=@RNA_EXON_ONLY/mRNA.biotype=@RNA_EXON_ONLY/exon
    # pseudogene SUB pseudogene/+pseudogenic_transcript/+exon
    # pseudogene/exon SUB pseudogene/+pseudogenic_transcript/exon
    # addd exon/add cds

    # gene	SUB	pseudogene/+pseudogenic_transcript/+exon	# MT related (ie MT ORI in flybase)
    # gene/cds	SUB	gene/+mRNA/exon	# MT related coding genes, should use both rules, perhaps use FIX rule?
    # gene/cds	SUB	gene/+mRNA/CDS

  def run_subsitutions(self, ctx, re_ctx):
    if ctx is None:
      return None
    if self._action is None:
      return None
    #
    it = ctx
    for ait in reversed(self._action):
      aa = ait["action"]
      if aa == "add":
        continue
      if aa == "rename":
        self.update_node(it, ait, re_ctx)
      #do nothing for exclude and copy_leaf
      it = it.get("_PARENTCTX")
    return

  def update_node(self, node, action, re_ctx):
    if not node or not action:
      return
    _type = self.from_rectx(action["type"], re_ctx)
    node["_TYPE"] = _type
    used_quals = node.get("_RULESDATA")["_ALL"].get("USEDQUALS")
    if used_quals is not None:
      for q, v in action.get("quals",{}).items():
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
    return

  def from_rectx(self, x, re_ctx=None):
    if not re_ctx or x is None: return x
    if type(x) != str or len(x) < 2: return x
    if x[0] == "@":
      x = x[1:]
      x = re_ctx.get(x, x)
    return x


