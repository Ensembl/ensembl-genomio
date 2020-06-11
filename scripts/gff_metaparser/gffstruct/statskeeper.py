import json
import sys

from collections import defaultdict, OrderedDict

from .basekeeper import BaseKeeper

class StatsKeeper(BaseKeeper):
  # storing stats
  def __init__(self):
    self._data = defaultdict(lambda: defaultdict(OrderedDict))

  def add(self, rule_name, context):
    print("adding stats for rule %s context %s" % (rule_name, context), file = sys.stderr)
    pass

  def summary(self, rule_name):
    pass

  def dump(self, out_file, detailed = None):
    pass

