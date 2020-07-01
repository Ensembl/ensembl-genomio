import json

from collections import defaultdict, OrderedDict

from .basekeeper import BaseKeeper

class FannKeeper(BaseKeeper):
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
      if not force:
        if not self.list_has(top, value):
          top.append(value)
      else:
        top.clear()
        top.append(value)
      return
    for k, v in value.items():
      if not v:
        continue
      if k not in top:
        top[k] = v
        continue
      if type(v) == list:
        if not force:
          top[k] += v
          top[k] = self.uniq_list(top[k])
        else:
          top[k] = [ v ]
        continue
      if not force:
        continue
      top[k] = v
    return

  def uniq_list(self, lst):
    # slow
    if len(lst) <= 1:
      return lst
    return [ json.loads(p) for p in frozenset(map(lambda x: json.dumps(x, sort_keys = True), lst)) ]

  def list_has(self, lst, obj):
    # slow
    if not lst or obj is None or len(lst) < 1:
      return False
    prj = json.dumps(obj, sort_keys = True)
    return prj in frozenset(map(lambda x: json.dumps(x, sort_keys = True), lst))


  def dump_json(self, out_file, maps = None, dump_filter=None):
    if not out_file:
      return
    vals = []
    for tag in sorted(self._data.keys()):
      vals += list(filter(dump_filter, self._data[tag].values()))
    json.dump(vals, out_file, indent = 2, sort_keys = True)

  def dump(self, out_file, maps = None, dump_filter=None):
    self.dump_json(out_file, maps, dump_filter)
