import copy
import sys

from .valid import ValidRule

class FixAction:
  def __init__(self, raw, rule, always_generate_new_ids = False):
     self._raw = raw
     self._rule_name = rule.NAME
     self._rule_lineno = rule._lineno
     self._always_gen_id = always_generate_new_ids
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
    possible_leaf = raw_splitted[-1]
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
        if id(p) != id(possible_leaf):
          self._action = False
          print("copy_leaf(!) used not with leaf in %s" % (self), file=sys.stderr)
          return None
        p = p[1:]
        self._copy_leaves += 1
      _type, *quals = p.replace(":",",").split(",")
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
    """
      copy node if updating parentctx:
        if adding or changing parentctx (skipping)
        remove/set _ISLEAF when needed
      return dict added nodes {id(node):id}
      mark removed/updated leaves as {id(node):none}
      mark leaves as to be removed from old leaves
    """
    if ctx is None:
      return None
    if self._action is None:
      return None
    depth = ctx.get("_DEPTH")
    if depth != len(self._action) - self._additions:
      fulltag = ctg.get("_FULLTAG")
      print("unbalanced number of tag %s and actions %s. skipping", fulltag, self, file = sys.stderr)
      return None

    re_ctx = ctx.get("_RECTX") and ctx["_RECTX"].groupdict() or None
    self.run_subsitutions(ctx, re_ctx)

    # data to use by "add"
    adepth = len(self._action or []) # action depth to use by "add"

    # add / remove, updating parent ctx
    new_nodes = {}
    prev, it = None, ctx
    for ait in reversed(self._action):
      is_leaf, parent = False, None
      if it:
        is_leaf = it.get("_ISLEAF", False)
        parent = it.get("_PARENTCTX")
      aa = ait["action"]
      if aa == "exclude":
        if is_leaf:
          self.del_node_if_leaf(it, new_nodes)
          if parent:
            parent = self.copy_node(parent, new_nodes)
            parent["_ISLEAF"] = True
        elif prev:
          prev = self.copy_node(prev, new_nodes)
          prev["_PARENTCTX"] = parent
        #
        prev, it, adepth = prev, parent, adepth - 1
        continue
      elif aa == "add":
        # x x N x x
        insert_after = bool(it)
        if insert_after:
          # copy (for the first time)
          it = self.copy_node(it, new_nodes, clean = False)
          _src_id = self.gen_and_update_id(it, parent = parent, depth = adepth)
          # copy to add
          _src = it
          _it = self.copy_node(_src, new_nodes, clean = True, force = True)
          _it["_PARENTCTX"] = _src
          # leaf or update kid
          if prev is not None:
            prev = self.copy_node(prev, new_nodes)
            prev["_PARENTCTX"] = _it
            _it["_ISLEAF"] = False
          else:
            _it["_ISLEAF"] = True
          # update node
          _gen_id = not _it.get("_ISLEAF", False) or self._always_gen_id
          _adata = self.copy_action(ait, gen_id = _gen_id, src_id = _src_id, depth = adepth)
          self.update_node(_it, _adata)
          #
          prev, it, adepth = _it, it, adepth - 1
          continue
        else: # insert before
          # force modify prev, otherwise copy the whole chain
          _src_id = self.gen_and_update_id(prev, parent = parent, depth = adepth)
          _src = prev
          # copy to add
          _it = self.copy_node(_src, new_nodes, clean = True, force = True)
          _it["_PARENTCTX"] = None
          _it["_ISLEAF"] = False
          # force update prev, otherwise, copy the whole chain
          prev["_PARENTCTX"] = _it
          #
          _adata = self.copy_action(ait, gen_id = True, src_id = _src_id, depth=adepth)
          self.update_node(_it, _adata)
          #
          prev, it, adepth = _it, it, adepth - 1
          continue
      elif aa == "copy_leaf":
        if is_leaf:
          it = self.copy_node(it, new_nodes, keep_leaf = True, clean = True)
          self.update_node(it, ait, re_ctx)
      #
      prev, it, adepth = it, parent, adepth - 1
    #
    return new_nodes

  def copy_action(self, action, gen_id = False, src_id = None, depth = 0, type = None):
    data = action.copy()
    data["quals"] = data.get("quals", {}).copy()
    if type: action["type"] = type
    if gen_id:
      data["quals"].update({ "ID" : self.new_id(src_id = src_id, depth = depth) })
    return data

  def new_id(self, src_id = None, depth = 0):
    if depth < 0: depth  = 100 - depth
    return "%s_dfd_%s" % (src_id, depth)

  def copy_node(self, node, new_nodes, keep_leaf = False, clean = False, force = False):
    if not force and node.get("_ISCOPY"):
      return node
    ncopy = node.copy()
    old_rules_data = ncopy.get("_RULESDATA")
    if old_rules_data is not None:
      ncopy["_RULESDATA"] = copy.deepcopy(old_rules_data)
    ncopy["_ISCOPY"] = True
    new_nodes[id(ncopy)] = ncopy
    # clean
    if clean:
      ncopy.get("_RULESDATA")["_ALL"]["USEDQUALS"] = {}
    # mark old leaf as unused
    if not keep_leaf:
      self.del_node_if_leaf(node, new_nodes)
    return ncopy

  def del_node_if_leaf(self, node, new_nodes):
    if node.get("_ISLEAF", False):
      new_nodes[id(node)] = None

  def run_subsitutions(self, ctx, re_ctx=None):
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

  def update_node(self, node, action, re_ctx = None):
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

  @classmethod
  def update_id(cls, node, id):
    if not node:
      return
    node["_ID"] = id
    used_quals = node.get("_RULESDATA")["_ALL"].get("USEDQUALS")
    if used_quals is None:
      node.get("_RULESDATA")["_ALL"]["USEDQUALS"] = {}
      used_quals = node.get("_RULESDATA")["_ALL"].get("USEDQUALS")
    used_quals.update({"id":("ID", id)})

  def gen_and_update_id(self, it, parent = None, depth = 0):
    _id = it.get("_ID")
    if _id:
      return _id
    if depth < 0: depth = 100 - depth
    _pid = parent and parent.get("_ID")
    if _pid:
      _id = "%s_dfcd_%s" % (_pid, depth) # derived feature from child at depth
    else:
      # use locus tag
      _id = "%s:%s_%s%s_dfd_%s" % (it["_SEQID"], it["_START"], it["_END"], it["_STRAND"] != "1" and it["_STRAND"] or "", depth)
    self.update_id(it, _id)

  def from_rectx(self, x, re_ctx=None):
    if not re_ctx or x is None: return x
    if type(x) != str or len(x) < 2: return x
    if x[0] == "@":
      x = x[1:]
      x = re_ctx.get(x, x)
    return x

