# valid structures

import sys

from .basestruct import *
from .rules import *

# VALID STRUCTURES class


class MetaParserStructures(BaseStructures):
  # ugly class names, rename to something more meaningful
  KNOWN_RULES = [
    AliasRule,
    GffRule,
    GffSubRule,
    GffForceSubRule,

    JsonRule,
    JsonSubRule,
    JsonForceSubRule,
    JsonAppendRule,
    JsonIdRule, # ? do we need it
  ]

  class Structure:
    def __init__(self, tag = ""):
      self.tag = tag

