from .base import BaseRule

class JsonRule(BaseRule):
  NAME = "JSON"
  _RULES = BaseRule.RulesType()

class JsonSubRule(JsonRule):
  NAME = "JSON_SUB"
  _RULES = BaseRule.RulesType()

class JsonForceSubRule(JsonSubRule):
  NAME = "JSON_FORCE_SUB"
  _RULES = BaseRule.RulesType()

class JsonAppendRule(JsonRule):
  # only one '@' can be used in the output path to mark array to  append to
  NAME = "JSON_APPEND"
  _RULES = BaseRule.RulesType()

class JsonIdRule(JsonRule):
  # do we need it at all
  NAME = "JSON_ID"
  _RULES = BaseRule.RulesType()


