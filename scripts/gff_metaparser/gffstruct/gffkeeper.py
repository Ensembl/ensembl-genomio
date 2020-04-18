from collections import defaultdict

class GffKeeper:
  # storing/building gff structures
  def __init__(self):
    self._data = defaultdict(lambda: defaultdict(dict))



