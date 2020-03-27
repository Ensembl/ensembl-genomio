from .base import BaseRule

class GffRule(BaseRule):
  NAME = "GFF"
  _RULES = BaseRule.RulesType()

class GffSubRule(GffRule):
  NAME = "GFF_SUB"
  _RULES = BaseRule.RulesType()

class GffForceSubRule(GffSubRule):
  NAME = "GFF_FORCE_SUB"
  _RULES = BaseRule.RulesType()

