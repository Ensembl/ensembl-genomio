import sys

from collections import defaultdict

from .base import BaseRule
from .valid import ValidRule
from .fix_action import FixAction


class SubRule(ValidRule):
  NAME = "SUB"
  _RULES = BaseRule.RulesType()

  def prepare_actions(self):
    self._actions = None
    raw = [ x.strip() for x in " ".join(self._actions_raw).split() if x.strip() ]
    self._actions = list(filter(lambda x: bool(x), [ FixAction(r, self) for r in raw ]))

  def process(self, context, re_context = None):
    # super(ValidRule) to fill stats and context
    ValidRule.process(context, re_context, name_override = self.NAME)
    # cations part
    actions = context.get("_RULESDATA")[self.NAME].get("ACTIONS")
    if actions is None:
      context.get("_RULESDATA")[self.NAME]["ACTIONS"] = list()
      actions = context.get("_RULESDATA")[self.NAME].get("ACTIONS")
    actions += self._actions

  @classmethod
  def run_postponed(cls, context, name_override = None):
    name_to_check = name_override or cls.NAME
    # super(ValidRule) to fill initial data
    ValidRule.run_postponed(context, name_override = name_to_check)

    # shared between actions and chains / leaves
    nodes_data = defaultdict(dict)
    #   new_nodes {id(node) : node}
    #
    for ctx in context.used_leaves():
      fulltag = ctx.get("_FULLTAG")
      check_quals = ctx.get("_RULESDATA")[name_to_check].get("USEDQUALS")
      if check_quals is None:
        continue
      actions = ctx.get("_RULESDATA")[cls.NAME].get("ACTIONS")
      if not actions:
        continue
      actions = frozenset(actions)
      #
      types_num = cls.action_types_num(actions)
      if types_num > 1:
        print("too many different action types (%s) for %s. skipping" %
               (types_num, fulltag), file = sys.stderr)
        continue
      #
      done = set()
      for a in actions:
        if id(a) in done: continue
        res = a.act(ctx, nodes_data)
        nodes_data["new_nodes"].update(res or {})
        done.add(id(a))
      #
    # alterations
    if not nodes_data["new_nodes"]:
      return
    new_nodes = nodes_data["new_nodes"]
    # drop leaves previously marked for delition
    ctx_leaves = context._useful_leaves
    print("\ncleaning old leaves", file=sys.stderr)
    for i in range(len(ctx_leaves))[::-1]:
      leaf = ctx_leaves[i]
      lid = id(leaf)
      if lid in new_nodes:
        if new_nodes[lid] is None or new_nodes[lid].get("_ISDELETED"):
          print("removin old leaf", lid, file=sys.stderr)
          ctx_leaves.pop(i)

    # update ctx and leaves
    for _id, node in new_nodes.items():
      if node is None:
        continue
      if node.get("_ISDELETED"):
        continue
      context.prev.append(node)
      if node.get("_ISLEAF"):
        print("adding leaf node with id", _id, file=sys.stderr)
        ctx_leaves.append(node)
    # find ids that met once
    seen_ids = []
    for leaf in ctx_leaves:
      seen_ids += context.get_to_root(node = leaf,
        getter = lambda x: (FixAction.n_id(x), id(x), x.get("_NOIDUPDATE")) ) or []
    seen_ids = set(seen_ids)
    print("seen_ids", seen_ids, file=sys.stderr)
    seen_counts = defaultdict(int)
    for _id, _ptr, _flag in seen_ids:
      if _id:
        seen_counts[_id] +=1
    met_once = frozenset([_id for _id, _c in seen_counts.items() if _c == 1])
    print("met_once", met_once, file=sys.stderr)
    #met_once = frozenset()
    print("met_once", met_once, file=sys.stderr)
    # iterate, check if _NOIDUPDATE is not none and false
    seen = defaultdict(int)
    seen.update({_id:1 for _id, _ in seen_counts.items()})
    print("seen_ids", seen, file=sys.stderr)
    updated_obj_ids = set()
    for leaf in ctx_leaves:
      context.run_to_root(node = leaf,
        updater = lambda x: cls.id_updater(x, seen, met_once, updated_obj_ids)
      )
    print("seen_ids", seen, file=sys.stderr)
    print("updated obj ids", updated_obj_ids, file=sys.stderr)
    return

  @classmethod
  def id_updater(cls, node, seen_ids, met_once, updated_obj_ids):
    if "_NOIDUPDATE" not in node or node["_NOIDUPDATE"]:
      return
    _id = FixAction.n_id(node)
    if not _id or _id in met_once:
      return
    cls.update_seen_id(node, seen_ids, updated_obj_ids)

  @classmethod
  def update_seen_id(cls, node, seen_ids, updated):
    if id(node) in updated:
      return
    print("updating", id(node), file=sys.stderr)
    updated.add(id(node))
    _id = FixAction.n_id(node)
    if not _id:
      return None
    seen_ids[_id] += 1
    if seen_ids[_id] == 1:
      return
    _id = "%s_%s" % (_id, seen_ids[_id])
    FixAction.update_id(node, _id)

  @classmethod
  def action_types_num(cls, actions):
    _w_ex = len(list(filter(lambda a: a._exclusions > 0, actions)))
    _w_add = len(list(filter(lambda a: a._additions > 0, actions)))
    # ignore copy leaves here
    actions_types_num = sum([
        _w_ex > 0,
        _w_add > 0,
        (len(actions) - _w_ex - _w_add) > 0,
    ])
    return actions_types_num


class SpellRule(BaseRule):
  NAME = "SPELL"
  _RULES = BaseRule.RulesType()

