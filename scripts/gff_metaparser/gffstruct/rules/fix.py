import sys

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

    new_nodes = dict() # {id(node) : node}
    for ctx in context.used_leaves():
      check_quals = ctx.get("_RULESDATA")[name_to_check].get("USEDQUALS")
      if check_quals is None: continue
      #
      actions = ctx.get("_RULESDATA")[cls.NAME].get("ACTIONS")
      if not actions: continue
      #
      _w_ex = len(list(filter(lambda a: a._exclusions > 0, actions)))
      _w_add = len(list(filter(lambda a: a._additions > 0, actions)))
      actions_types_num = sum([
          _w_ex > 0,
          _w_add > 0,
          (len(actions) - _w_ex - _w_add) > 0,
      ])
      if actions_types_num > 1:
        fulltag = ctx.get("_FULLTAG")
        print("too many action types for %s. skipping" % fulltag, file = sys.stderr)
        continue
      #
      done = set()
      for a in actions:
        if id(a) in done: continue
        done.add(id(a))
        res = a.act(ctx)
        new_nodes.update(res or {})

    # alterations
    if not new_nodes: return
    # drop unused leaves
    ctx_leaves = context._useful_leaves
    for i in range(len(ctx_leaves))[::-1]:
      leaf = ctx_leaves[i]
      if id(leaf) in new_nodes and new_nodes[id(leaf)] is None:
        ctx_leaves.pop(i)
    # append new nodes
    for _id, node in new_nodes.items():
      if node is None: continue
      context.prev.append(node)
      if node.get("_ISLEAF"): ctx_leaves.append(node)
    #
    return


class SpellRule(BaseRule):
  NAME = "SPELL"
  _RULES = BaseRule.RulesType()

