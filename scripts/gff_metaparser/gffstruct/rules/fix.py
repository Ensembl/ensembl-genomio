import sys

from .base import BaseRule
from .valid import ValidRule

class SubRule(ValidRule):
  NAME = "SUB"
  _RULES = BaseRule.RulesType()

  @classmethod
  def run_postponed(cls, context, name_override = None):
    name_to_check = name_override or cls.NAME
    # super(ValidRule) to fill initial data
    ValidRule.run_postponed(context, name_override = name_to_check)

    for ctx in context.used_leaves():
      check_quals = ctx.get("_RULESDATA")[name_to_check].get("USEDQUALS")
      if check_quals is None:
        continue
      re_ctx = ctx.get("_RECTX")
      if re_ctx:
        re_ctx = re_ctx.groupdict() or None
        print("run_postponed for sub rectx:", re_ctx, ctx.get("_FULLTAG"), file=sys.stderr)


class SpellRule(BaseRule):
  NAME = "SPELL"
  _RULES = BaseRule.RulesType()

