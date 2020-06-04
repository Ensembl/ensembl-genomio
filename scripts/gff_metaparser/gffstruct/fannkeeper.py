import json

from collections import defaultdict, OrderedDict

class FannKeeper:
  # storing result functional annotation object
  def __init__(self):
    self._data = defaultdict(lambda: defaultdict(OrderedDict))

  def add(self, obj_tag, obj_id, path, value, force = False):
    if not value:
      return

    # find a place to insert to
    top = self._data[obj_tag][obj_id]
    if path:
      for p in path.split("/"):
        if type(p) == list:
          Except("To many list levels for %s:%s/%s" %(obj_tag, obj_id, path))
        p_is_list = False
        if p.startswith("@"):
          p_is_list = True
          p = p[1:]
        if p not in top:
          if p_is_list:
            top[p] = list()
          else:
            top[p] = dict()
        top = top[p]

    # update top
    if type(top) == list:
      top.append(value)
      return
    for k, v in value.items():
      if not v:
        continue
      if k not in top:
        top[k] = v
        continue
      if type(v) == list:
        top[k] += v
        top[k] = self.uniq_list(top[k])
        continue
      if not force:
        continue
      top[k] = v
    return

  def uniq_list(self, lst):
    # slow
    if len(lst) <= 1:
      return lst
    return [ json.loads(x) for x in frozenset(map(lambda x: json.dumps(x, sort_keys = True), lst)) ]

  def dump_json(self, out_file, maps = None, dump_filter=None):
    if not out_file:
      return
    vals = []
    for tag in sorted(self._data.keys()):
      vals += list(filter(dump_filter, self._data[tag].values()))
    json.dump(vals, out_file, indent = 2, sort_keys = True)
