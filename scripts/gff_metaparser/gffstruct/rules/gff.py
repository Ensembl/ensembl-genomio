from .base import BaseRule

class GffRule(BaseRule):
  # IDs can be substituted, but can be  ommited only for leaves
  # parents can't be ommited from GFF, paresnt, phase and strand  will be preserved in GFF
  # though parents should be excluded from qulifiers, they'll be added automatically based on parent feature ID
  NAME = "GFF"
  _RULES = BaseRule.RulesType()

class GffSubRule(GffRule):
  NAME = "GFF_SUB"
  _RULES = BaseRule.RulesType()

class GffForceSubRule(GffSubRule):
  NAME = "GFF_FORCE_SUB"
  _RULES = BaseRule.RulesType()

