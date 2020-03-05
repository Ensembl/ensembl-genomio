from .base import BaseRule

class UnseenRule(BaseRule):
  @classmethod
  def process(cls, struct, noconfig=True):
    print("no matching pattern for %s" % struct.tag, file=sys.stderr)
    pass

