# rules classes
from collections import defaultdict

class BaseRule:
  NAME = ""
  _const_patterns = {}
  _regex_patterns = []

  def __init__(self, pattern, actions, lineno):
    pass

  @classmethod
  def const_patterns(cls):
    return cls._const_patterns.keys()

  @classmethod
  def regex_patterns(cls):
    return list(map(lambda x: x[0], cls._regex_patterns))

  def const_match(self, tag):
    pass

  def regex_match(self, tag):
    pass

  def process(self, tag, re_context = None):
    pass

  def dump(self):
    pass

