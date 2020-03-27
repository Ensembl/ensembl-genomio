# valid structures

import sys

from .basestruct import *
from .rules import *

# VALID STRUCTURES class

class ValidStructures(BaseStructures):
  KNOWN_RULES = [
    AliasRule,
    GffRule,
    GffSubRule,
    GffForceSubRule,

    JsonRule,
    JsonSubRule,
    JsonForceSubRule,
    JsonAppendRule,
    JsonIdRule,
  ]

  class Structure:
    def __init__(self, tag = ""):
      self.tag = tag

