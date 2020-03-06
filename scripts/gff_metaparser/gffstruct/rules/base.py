# rules classes

import sys
from collections import defaultdict
from copy import deepcopy

class BaseRule:
  NAME = ""
  _const_patterns = {}
  _regex_patterns = []
  _RULES = dict()

  def __init__(self, pattern, actions, lineno):
    self._pattern = pattern
    self._actions = actions
    self._lineno = lineno
    #BaseRule.update_rules(self)
    self.update_rules()

  def update_rules(self):
   pat = self._pattern.strip().lower()
   if pat in self._RULES:
     print("already seen pattern %s at %d for %s at %d" % (
         self._pattern, self._lineno, self.NAME, self._RULES[pat]._lineno
       ), file=sys.stderr)
   else:
     self._RULES[pat] = self

  @classmethod
  def const_patterns(cls):
    return deepcopy(cls._const_patterns)

  @classmethod
  def regex_patterns(cls):
    return deepcopy(cls._regex_patterns)

  def const_match(self, tag):
    pass

  def regex_match(self, tag):
    pass


  def process(self, tag, re_context = None):
    pass

  def dump(self):
    pass

