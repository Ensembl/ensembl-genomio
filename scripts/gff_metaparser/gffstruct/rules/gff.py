from .base import BaseRule

class GffRule(BaseRule):
  # IDs can be substituted, but can be  ommited only for leaves
  # parents can't be ommited from GFF, paresnt, phase and strand  will be preserved in GFF
  # though parents should be excluded from qulifiers, they'll be added automatically based on parent feature ID
  NAME = "GFF"
  _RULES = BaseRule.RulesType()

  def prepare_actions(self):
     pass

  def process(self, context, re_context = None):
    print("processing %s for %s with match groups %s" % (
        context.tag(), self.NAME, str(re_context and re_context.groupdict() or None)
      ))
    pass


class GffSubRule(GffRule):
  NAME = "GFF_SUB"
  _RULES = BaseRule.RulesType()

