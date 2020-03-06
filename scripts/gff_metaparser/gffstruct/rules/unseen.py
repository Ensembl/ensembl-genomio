import sys

from .base import BaseRule

class UnseenRule(BaseRule):
  NAME = "UNSEEN"
  _RULES = BaseRule.RulesType()

  @classmethod
  def process(cls, struct, noconfig=True):
    print("no matching pattern for %s" % struct.tag, file=sys.stderr)
    pass

