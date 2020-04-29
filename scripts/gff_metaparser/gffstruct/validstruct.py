# valid structures

import sys

from .basestruct import *
from .rules import *

# VALID STRUCTURES class

class ValidStructures(BaseStructures):
  KNOWN_RULES = [
    AliasRule,
    IgnoreRule,
    ValidRule,
    ValidIfRule,
    FixRule,
    ForceFixRule,
    SubRule,
    ForceSubRule,
    SetRule,
  ]

