import sys

from .base import BaseRule

class UnseenRule(BaseRule):
  NAME = "UNSEEN"
  _RULES = BaseRule.RulesType()

  @classmethod
  def process(cls, context, re_context = None):
    # add stats
    context.global_context.add(cls.NAME, context)

