from .base import BaseRule

class ValidRule(BaseRule):
  NAME = "VALID"
  _RULES = BaseRule.RulesType()

  @classmethod
  def process(cls, context, re_context = None):
    # add stats
    context.global_context.add(cls.NAME, context)


class ValidIfRule(ValidRule):
  NAME = "VALID_IF"
  _RULES = BaseRule.RulesType()
  # store gene.id/mrna.id/_feature at global context for checking

  @classmethod
  def prepare_postponed(cls, context):
    pass

  @classmethod
  def run_postponed(clsf, context):
    pass
