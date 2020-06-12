import json
import sys

from collections import defaultdict, OrderedDict

from .basekeeper import BaseKeeper
from .prefixtree import PfxTr

class StatsKeeper(BaseKeeper):
  # storing stats
  def __init__(self, detailed = None, no_longest_pfx = False):
    self._data = defaultdict(lambda: defaultdict(dict))
    self._detailed = detailed
    self._no_longest_pfx = no_longest_pfx

  def add(self, rule_name, context):
    _fulltag = context.get("_FULLTAG")
    tagstat = self._data[rule_name][_fulltag]

    if "counts" not in tagstat:
      tagstat["counts"] = 0
    tagstat["counts"] += 1

    type_id_coords = [ self.type_id_coords_from_ctx(x) for x in context.prev + [context] ]

    if self._detailed:
      print("rule %s for %s: %s " % (rule_name, _fulltag, type_id_coords), file = self._detailed)

    if not self._no_longest_pfx:
      if "pfx" not in tagstat:
        tagstat["pfx"] = defaultdict(PfxTr)
      for _type, _id, _ in type_id_coords:
        if _id and _id != ".":
          tagstat["pfx"][_type].add(_id)

  def summary(self, rule_name, out_file = None):
    if rule_name not in self._data:
      return None
    out = []
    for _tag in self._data[rule_name]:
      tagstat = self._data[rule_name][_tag]
      counts = str(tagstat.get("counts"))
      chain = str([ { tp : pfx.get_max() } for tp, pfx in tagstat.get("pfx", {}).items() ])
      out_line = "\t".join(["#stats", rule_name, _tag, counts, chain])
      if out_file:
        print(out_line, file = out_file)
      else:
        out.append(out_line)
    if out_file:
      return True
    return "\n".join(out)

  def dump(self, out_file):
    for _rule in self._data:
      self.summary(_rule, out_file)

  def type_id_coords_from_ctx(self, context):
    if not context:
      return ".", ".", "."
    _type, _id, _seqid, _start, _end, _strand = map(lambda x: context.get("_"+x) or ".", "TYPE ID SEQID START END STRAND".split())
    _strand = _strand == "." and "." or int(_strand) < 0 and "-" or "+"
    _coords = "%s:%s:%s-%s" % (_seqid, _strand, _start, _end)
    return _type, _id, _coords
