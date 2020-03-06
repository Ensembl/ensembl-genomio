# rules classes

import sys

from collections import defaultdict
from copy import deepcopy


def _RulesType():
  return {
    "const_match" : defaultdict(list),
    "_regex_match_pre" : defaultdict(list),
    "regex_match" : list(),
  }


class BaseRule:
  NAME = ""
  _RULES = _RulesType()

  @classmethod
  def RulesType(cls):
    return _RulesType()

  @classmethod
  def AliasSymbol(cls):
    return r'@'

  def __init__(self, pattern, actions, lineno):
    self._pattern = pattern
    self._actions = actions
    self._lineno = lineno
    self.update_rules()

  def update_rules(self):
    pat = self._pattern.strip().lower()
    rules = self._RULES["const_match"]
    if BaseRule.AliasSymbol() in pat:
      rules = self._RULES["_regex_match_pre"]

    if pat in rules:
      print("already seen pattern %s at lines %s for %s at %d" % (
          self._pattern, ",".join(map(lambda r: str(r._lineno), rules[pat])),
          self.NAME, self._lineno,
        ), file=sys.stderr)
    rules[pat].append(self)

  @classmethod
  def mature_regex(cls, aliases_cls):
    if not aliases_cls:
      return
    for pat, rules in cls._RULES["_regex_match_pre"].items():
      for rule in rules:
        raw_pat = rule._pattern.strip()
        re_pat, re_compiled = aliases_cls.mature_regex(raw_pat)
        if re_pat == None:
          continue
        if re_pat == raw_pat:
          cls._RULES["const_match"][pat].append(rule)
          print("cannot mature pattern %s (raw: %s) for %s (line %d)" % (
               pat, raw_pat, cls.NAME, rule._lineno,
            ), file=sys.stderr)
        else:
            cls._RULES["regex_match"].append((re_pat, re_compiled, rule))

  @classmethod
  def const_patterns(cls):
    return deepcopy(cls._RULES["const_match"])

  @classmethod
  def regex_patterns(cls):
    return deepcopy(cls._RULES["regex_match"])

  def process(self, tag, re_context = None):
    pass

