import copy

class WalkContext:
  def __init__(self, tag = "", gff_keeper = None, global_context = None, ctg_len_inferer = None):
    self.data = dict()
    self._tag = tag
    self.gff_keeper = gff_keeper
    self.global_context = global_context
    self.ctg_len = ctg_len_inferer
    self.processed_rules = []
    self.prev = []
    self._top = None
    self.fixes = []

  def snap(self):
    # shallow data copy
    self.prev.append({"tag" : self._tag, "data" : copy.copy(self.data)})
    return self.prev[-1]

  def top(self, *feature):
    if len(feature) == 0:
      return self._top
    if feature[0]:
      self._top = feature[0]

  def tag(self, *val):
    if len(val) == 0:
      return self._tag
    self._tag = str(val[0])

  def update(self, *key_val, force_clean=False, **kwargs):
    # can have either dict as key or string with non-empty val
    if len(key_val) == 1 and isinstance(key_val[0], dict):
      for k, v in key_val[0].items():
        self.update(k,v)
    elif len(key_val) == 2:
      key, val = key_val
      if val is not None:
        self.data[key] = val
      elif force_clean and k in self.data:
        del self.data[k]
    # update from **kwargs
    for k, v in kwargs.items():
      self.update(k,v)

  def get(self, key, default = None):
    # global ???
    if key not in self.data:
      return default
    return self.data[key]

  def __getitem__(self, key):
    return self.get(key)

  def update_processed_rules(self, lst):
    self.processed_rules += lst

  def add_fix(self, fix):
    pass

  def merge_fixes(self):
    pass

