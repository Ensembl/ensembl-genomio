from .base import BaseRule

class SubRule(BaseRule):
  NAME = "SUB"
  _RULES = BaseRule.RulesType()

class SpellRule(BaseRule):
  NAME = "SPELL"
  _RULES = BaseRule.RulesType()

