from .base import BaseRule

class SubRule(BaseRule):
  NAME = "SUB"
  _RULES = BaseRule.RulesType()

  @classmethod
  def process(cls, context, re_context = None):
    # add stats
    context.global_context.add(cls.NAME, context)


class SpellRule(BaseRule):
  NAME = "SPELL"
  _RULES = BaseRule.RulesType()

