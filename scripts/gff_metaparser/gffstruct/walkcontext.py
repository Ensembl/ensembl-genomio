import copy

class WalkContext:
  def __init__(self, tag = "", global_context = None, ctg_len_inferer = None):
    self.data = dict()
    self.tag = tag
    self._gctx = global_context
    self._ctg_len = ctg_len_inferer
    self.processed_rules = []
    self.prev = []
    self.top = None
    self.fixes = []

  def snap():
    # shallow data copy
    self.prev.append({"tag" : tag, "data" : copy.copy(self.data)})

  def top(*feature):
    if len(feature) == 0:
      return self.top
    if feature[0]:
      self.top = feature[0]

  def tag(*val):
    if len(val) == 0:
      return self.tag
    self.tag = str(val[0])

  def update(self, *key_val, **kwargs):
    # can have either dict as key or string with non-empty val
    if isinstance(key_val, dict):
      for k, v in key_val.items():
        self.update(k,v)
    elif len(key_val) == 2:
      key, val = key_val
      if val is not None:
        data[key] = val
    # update from **kwargs
    for k, v in kwargs.items():
      self.update(k,v)

  def get(self, key):
    # global ???
    if key not in data:
      return None
    return data[key]
    pass

  def __getitem__(self, key):
    return self.get(key)

  def update_processed_rules(self, lst):
    self.processed_rules += lst

  def add_fix(self, fix):
    pass

  def merge_fixes(self):
    pass

