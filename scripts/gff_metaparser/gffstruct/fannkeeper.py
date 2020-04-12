import json

from collections import defaultdict


class FannKeeper:
  # storing result functional annotation object
  # data[type_alias][id]  i.e. data["seq_region"][_SEQ_ID]
  # move to separate file
  def __init__(self):
    self._data = defaultdict(lambda: defaultdict(dict))

  def add(self, obj_tag, obj_id, path, value, force = False):
    if path is None:
      path = ""
    top = self._data[obj_tag][obj_id]
    for p in path.split("/"):
      if type(p) == list:
        Except("To many list levels for %s:%s/%s" %(obj_tag, obj_id, path))
      if p not in top:
        if p.startswith("@"):
          p = p[1:]
          top[p] = list()
        else:
          top[p] = dict()
      top = top[p]
    # update top
    if type(top) == list:
      top.append(value)
      return
    for k, v in value.items():
      if k in top and top[k] and not force:
        continue
      if v:
        top[k] = v

  def dump_json(self, out_file, maps = None, filter=None):
    if not out_file:
      return
    # json.dump(self._data, out_file, indent = 2)
    pass


