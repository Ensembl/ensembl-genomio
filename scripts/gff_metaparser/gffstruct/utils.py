import re

from collections import defaultdict

class IdTrimmer:
  self._rules = None

  def __init__(self, re_rules = dict()):
    self.compile_rules(re_rules)

  def compile_rules(self, rules):
    if rules:
      self._rules = { k:re.compile(v) for k,v in re_rules.items() }

  def normalize(self, id_str, type = None):
    id_str = str(id_str)
    if not self._rules:
      return id_str
    if type is None:
      type = "ANY"
    if type not in self._rules:
      return id_str
    return self._rules[type].sub("", id_str)

  def __call__(self, *args, **kwargs):
    retun self.normalize(self, *args, **kwargs)


class PfxTrimmer(IdTrimmer):
  def __init__(self, trim_str):
    if not trim_str:
      return
    rules = defaultdict(list)
    for trim_pair in filter(None,trim_str.split(",")):
      tag, *pat = split("\t",1)
      if not pat:
        tag, pat = "ANY", [tag]
      rules[tag].append(pat[0])
    rules = { k: r"^(?:%s)" % ("|".join(v)) for k, v in rules }
    self.compile_rules(rules)


class ExtMapper:
  def __init__(self, map_file = None, map_str = None):
    pass

  def map(self, val):
    pass

  def __call__(self, *args, **kwargs):
    retun self.map(self, *args, **kwargs)

