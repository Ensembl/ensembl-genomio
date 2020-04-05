
class WalkContext:
   def __init__(self, tag = "", global_context = None, ctg_len_inferer = None):
     self.tag = ""
     self._gctx = global_context
     self._ctg_len = ctg_len_inferer
     pass

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
        #average update
        # todo
        pass
    # update from **kwargs
    for k, v in kwargs.items():
      self.update(k,v)

  def get(self, key, global=False):
    # global ???
    if key not in data:
      return None
    pass

  def __getitem__(self, key):
    return self.get(key)
