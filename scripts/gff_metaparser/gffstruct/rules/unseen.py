import sys

from .base import BaseRule

class UnseenRule(BaseRule):
  NAME = "UNSEEN"
  _RULES = BaseRule.RulesType()

  @classmethod
  def process(cls, context, noconfig=True):
    print("no matching pattern for %s" % context.tag(), file=sys.stderr)
    pass

