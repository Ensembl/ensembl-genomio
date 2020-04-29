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

    JsonRule,
    JsonSubRule,
  ]

