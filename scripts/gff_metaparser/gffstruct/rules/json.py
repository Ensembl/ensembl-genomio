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
  NAME = "JSON_APPEND"
  _RULES = BaseRule.RulesType()

class JsonIdRule(JsonRule):
  NAME = "JSON_ID"
  _RULES = BaseRule.RulesType()


