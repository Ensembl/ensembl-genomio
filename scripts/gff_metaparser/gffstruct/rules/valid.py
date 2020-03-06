from .base import BaseRule

class ValidRule(BaseRule):
  NAME = "VALID"
  _RULES = BaseRule.RulesType()

class ValidIfRule(ValidRule):
  NAME = "VALID_IF"
  _RULES = BaseRule.RulesType()

