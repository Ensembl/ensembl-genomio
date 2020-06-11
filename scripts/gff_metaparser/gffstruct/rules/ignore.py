from .base import BaseRule

class IgnoreRule(BaseRule):
  NAME = "IGNORE"
  _RULES = BaseRule.RulesType()

  @classmethod
  def process(cls, context, re_context = None):
    # add stats
    context.global_context.add(cls.NAME, context)
