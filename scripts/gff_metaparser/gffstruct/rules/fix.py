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

    for ctx in context.used_leaves():
      check_quals = ctx.get("_RULESDATA")[name_to_check].get("USEDQUALS")
      if check_quals is None: continue

      actions = ctx.get("_RULESDATA")[cls.NAME].get("ACTIONS")
      if not actions: continue

      re_ctx = ctx.get("_RECTX") and ctx["_RECTX"].groupdict() or None
      print("run_postponed for sub rectx:", re_ctx, ctx.get("_FULLTAG"), actions, file=sys.stderr)
      # run_postponed for sub rectx: {'MRNA': 'mRNA', 'CDS': 'CDS'} gene/mRNA/CDS ['gene/transcript.biotype=@MRNA/@CDS']
    return


class SpellRule(BaseRule):
  NAME = "SPELL"
  _RULES = BaseRule.RulesType()

