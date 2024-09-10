# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import copy
import re
import sys


class FixAction:
    _split_camel_re = re.compile(r"([a-z])([A-Z])")

    def __init__(self, raw, rule, always_generate_new_ids=False):
        self._raw = raw
        self._rule_name = rule.NAME
        self._rule_lineno = rule._lineno
        self._always_gen_id = always_generate_new_ids
        #
        self._additions = 0
        self._exclusions = 0
        self._copy_leaves = 0
        #
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
        return "action: %s of %s at %s" % (self._raw, self._rule_name, self._rule_lineno)

    def parse(self, raw):
        out = []
        raw_splitted = raw.split("/")
        possible_leaf = raw_splitted[-1]
        for p in raw_splitted:
            if p.strip() == "" or p == "-":
                out.append({"action": "exclude"})
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
            _type, *quals = p.replace(":", ",").split(",")
            _type = _type.strip()
            if _type == "":
                print("empty type for %s" % (self), file=sys.stderr)
                return None
            quals = dict([(kv + "=").split("=")[0:2] for kv in quals if kv])
            out.append(
                {
                    "action": action,
                    "type": _type,
                    "quals": quals,
                    "depth": len(out) + 1,
                }
            )
        return out or None

    def act(self, ctx, nodes_data):
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
        fulltag = ctx.get("_FULLTAG")
        depth = ctx.get("_DEPTH")
        if depth != self.active_depth():
            print("unbalanced number of tag %s and actions %s. skipping" % (fulltag, self), file=sys.stderr)
            return None
        # no mixed add & remove rules yet
        if self.action_types_num() > 1:
            print("too many action types for %s and action %s. skipping" % (fulltag, self), file=sys.stderr)
            return None

        # run substitutions on original nodes
        re_ctx = ctx.get("_RECTX") and ctx["_RECTX"].groupdict() or None
        self.run_subsitutions(ctx, re_ctx)

        new_nodes = {}
        # remove
        if self._exclusions > 0:
            self.run_exclusions(ctx, new_nodes)
            return new_nodes

        # copy leaves first
        if self._copy_leaves > 0:
            self.run_copy_leaves(ctx, new_nodes, re_ctx)
        # add
        if self._additions > 0:
            new_leaves = list(new_nodes.values())
            for node in new_leaves or [ctx]:
                self.run_additions(node, re_ctx, new_nodes, nodes_data)
            return new_nodes
        #
        return new_nodes

    def run_copy_leaves(self, ctx, new_nodes, re_ctx):
        prev, it = None, ctx
        for ait in reversed(self._action):
            is_leaf = it.get("_ISLEAF", False)
            parent = it.get("_PARENTCTX")
            #
            if ait["action"] == "copy_leaf":
                if is_leaf:
                    # copy node, keep leaf, add link to the source group
                    it = self.copy_node(it, new_nodes, keep_leaf=True, clean=True)
                    self.update_node(it, ait, re_ctx)
                    if self._always_gen_id:
                        src_id = self.n_id(parent)
                        if not src_id:
                            if parent:
                                src_id = self.new_id(parent["_TYPE"], ait["type"])
                                self.update_id(parent, src_id)
                            else:
                                src_id = ait["type"]
                        new_id = self.new_id(self.n_id(parent), ait["type"])
                        self.update_id(it, new_id)
                        it["_NOIDUPDATE"] = False
                #
            # copy only leaf, so
            break
            #
            prev, it = it, parent
            continue
        return

    def node2add_after(self, oit, action, prev, oprev, re_ctx, nodes_data, new_nodes):
        prev_nodes = nodes_data["prev_nodes"]
        aa, atype, adepth = action["action"], action["type"], action["depth"]
        if aa == "copy_leaf":
            aa = "rename"
        oit_id = oit is None and "None" or str(id(oit))
        pn_key = "%s_%s_%s" % (oit_id, aa, adepth)
        if pn_key in prev_nodes:
            new_node = prev_nodes[pn_key]
            return new_node
        #
        modify = aa == "add"
        new_node = self.copy_node(oit or oprev, new_nodes, clean=modify, force=True)
        self.del_node_if_leaf(oit, new_nodes)
        if modify:
            self.update_node(new_node, action, re_ctx)
            # gen id
            src_id = self.n_id(oit or oprev)
            if not src_id:
                if oit:
                    src_id = self.new_id(oit["_TYPE"], atype)
                    self.update_id(oit, src_id)
                elif oprev:
                    src_id = self.new_id(oprev["_TYPE"], atype)
                    self.update_id(oprev, src_id)
                else:
                    src_id = atype
            new_id = self.new_id(src_id, atype)
            self.update_id(new_node, new_id)
            new_node["_NOIDUPDATE"] = False
        else:
            if new_node.get("_NOIDUPDATE") is None:
                new_node["_NOIDUPDATE"] = True
        #
        if prev:
            prev["_PARENTCTX"] = new_node
            new_node["_ISLEAF"] = False
        else:
            new_node["_ISLEAF"] = True
        #
        prev_nodes[pn_key] = new_node
        return new_node

    def run_additions(self, ctx, re_ctx, new_nodes, nodes_data):
        # are different length allowed?
        # original and current iterators
        oit, oprev = ctx, None
        prev = None
        for ait in reversed(self._action):
            oparent = oit and oit.get("_PARENTCTX") or None
            it = self.node2add_after(oit, ait, prev, oprev, re_ctx, nodes_data, new_nodes)
            if prev:
                prev["_PARENTCTX"] = it
            prev = it
            if ait["action"] != "add":
                oit, oprev = oparent, oit
        #
        return

    def copy_node(self, node, new_nodes, keep_leaf=False, clean=False, force=False):
        if not force and node.get("_ISCOPY"):
            return node
        ncopy = node.copy()
        ncopy["_ISDELETED"] = False
        ncopy["_ISCOPY"] = not force
        new_nodes[id(ncopy)] = ncopy
        #
        old_rules_data = ncopy.get("_RULESDATA")
        if old_rules_data is not None:
            ncopy["_RULESDATA"] = copy.deepcopy(old_rules_data)
        # clean
        if clean:
            ncopy.get("_RULESDATA")["_ALL"]["USEDQUALS"] = {}
        # mark old leaf as unused
        if not keep_leaf:
            self.del_node_if_leaf(node, new_nodes)
        return ncopy

    def copy_chain(self, chain_start, stop_at, new_nodes, keep_leaf=False):
        if chain_start is None:
            return None
        if stop_at is None:
            return None
        #
        it, prev = chain_start, None
        while it:
            cp = self.copy_node(it, new_nodes)
            if prev:
                prev["_PARENTCTX"] = cp
            else:
                cp["_ISLEAF"] = True
            if it == stop_at:
                it = cp
                break
            prev, it = cp, cp.get("_PARENTCTX")
        #
        if not it:
            return None
        return it

    def del_node_if_leaf(self, node, new_nodes):
        if node and node.get("_ISLEAF", False):
            node_id = id(node)
            if new_nodes.get(node_id):
                node["_ISDELETED"] = True
            else:
                new_nodes[node_id] = None

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
            # do nothing for exclude and copy_leaf
            it = it.get("_PARENTCTX")
        return

    def run_exclusions(self, ctx, new_nodes):
        prev, it = None, ctx
        chain_start = ctx
        for ait in reversed(self._action):
            is_leaf = it.get("_ISLEAF", False)
            parent = it.get("_PARENTCTX")
            #
            if ait["action"] == "exclude":
                if is_leaf:
                    self.del_node_if_leaf(it, new_nodes)
                    if parent:
                        parent = self.copy_node(parent, new_nodes)
                        parent["_ISLEAF"] = True
                        chain_start = parent
                elif prev:
                    # copy the whole chain
                    prev = self.copy_chain(chain_start, prev, new_nodes)
                    prev["_PARENTCTX"] = parent
                #
                prev, it = prev, parent
            else:
                prev, it = it, parent
        return

    def update_node(self, node, action, re_ctx=None):
        if not node or not action:
            return
        _type = self.from_rectx(action["type"], re_ctx)
        node["_TYPE"] = _type
        used_quals = node.get("_RULESDATA")["_ALL"].get("USEDQUALS")
        if used_quals is not None:
            for q, v in action.get("quals", {}).items():
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
                    used_quals.update({_q: (q, v)})
        return

    @classmethod
    def n_id(cls, node):
        _id = node.get("_RULESDATA")["_ALL"].get("USEDQUALS", {}).get("id", [None, None])[1]
        if _id and type(_id) == list:
            _id = _id[0]
        if not _id:
            _id = node.get("_ID")
        if not _id:
            return None
        return _id

    @classmethod
    def update_id(cls, node, _id):
        if not node:
            return
        node["_ID"] = _id
        used_quals = node.get("_RULESDATA")["_ALL"].get("USEDQUALS")
        if used_quals is None:
            node.get("_RULESDATA")["_ALL"]["USEDQUALS"] = {}
            used_quals = node.get("_RULESDATA")["_ALL"].get("USEDQUALS")
        used_quals.update({"id": ("ID", _id)})

    @classmethod
    def new_id(cls, src_id, type=None, depth=0):
        if depth < 0:
            depth = 100 - depth
        _type = "df"
        if type:
            splitted_name = cls._split_camel_re.sub(r"\1_\2", type).split("_")
            _type = "df_" + "".join([s[0] for s in splitted_name if s])
        return "%s_%s" % (src_id, _type.lower())
        # return "%s_%sd%s" % (src_id, _type.lower(), depth)

    def __gen_and_update_id(self, it, parent=None, depth=0, action=None):
        _id = it.get("_ID")
        if _id:
            return _id
        _src_id = parent and parent.get("_ID")
        if _src_id:
            _src_id = "%s_ch" % _src_id  # child derived
        else:
            # use locus tag
            _src_id = "%s:%s_%s%s" % (
                it["_SEQID"],
                it["_START"] + 1,
                it["_END"],
                it["_STRAND"] != "1" and it["_STRAND"] or "",
            )
        _type = action and action.get("type") or None
        _id = self.new_id(src_id=_src_id, depth=depth, type=_type)

        self.update_id(it, _id)

    def from_rectx(self, x, re_ctx=None):
        if not re_ctx or x is None:
            return x
        if type(x) != str or len(x) < 2:
            return x
        if x[0] == "@":
            x = x[1:]
            x = re_ctx.get(x, x)
        return x

    def active_depth(self):
        return len(self._action) - self._additions

    def action_types_num(self):
        actions = [self]
        _w_ex = len(list(filter(lambda a: a._exclusions > 0, actions)))
        _w_add = len(list(filter(lambda a: a._additions > 0, actions)))
        # ignore copy leaves here
        actions_types_num = sum(
            [
                _w_ex > 0,
                _w_add > 0,
                (len(actions) - _w_ex - _w_add) > 0,
            ]
        )
        return actions_types_num
