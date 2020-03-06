from .base import BaseRule

class FixRule(BaseRule):
  NAME = "FIX"
  _RULES = BaseRule.RulesType()

class ForceFixRule(FixRule):
  NAME = "FORCE_FIX"
  _RULES = BaseRule.RulesType()

class SubRule(FixRule):
  NAME = "SUB"
  _RULES = BaseRule.RulesType()

class ForceSubRule(SubRule):
  NAME = "FORCE_SUB"
  _RULES = BaseRule.RulesType()

class SetRule(FixRule):
  NAME = "SET"
  _RULES = BaseRule.RulesType()

