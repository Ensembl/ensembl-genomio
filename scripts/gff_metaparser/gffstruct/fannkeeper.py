import json

from collections import defaultdict


class FannKeeper:
  # storing result functional annotation object
  # data[type_alias][id]  i.e. data["seq_region"][_SEQ_ID]
  # move to separate file
  def __init__(self):
    self._data = defaultdict(dict)

  def add(self, path, value):
    pass

  def dump_json(self, out_file, filter=None):
    json.dump(json_out, out_file, indent = 2)
    pass


