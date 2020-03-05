from .base import BaseRule

class FixRule(BaseRule):
  NAME = "FIX"

class ForceFixRule(FixRule):
  NAME = "FORCE_FIX"

class SubRule(FixRule):
  NAME = "SUB"

class ForceSubRule(SubRule):
  NAME = "FORCE_SUB"

class SetRule(FixRule):
  NAME = "SET"

