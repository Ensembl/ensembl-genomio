from .base import BaseRule

class ValidRule(BaseRule):
  NAME = "VALID"
  _RULES = dict()

class ValidIfRule(ValidRule):
  NAME = "VALID_IF"
  _RULES = dict()

