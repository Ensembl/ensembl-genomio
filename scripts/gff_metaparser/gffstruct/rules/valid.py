from .base import BaseRule

import sys
from collections import defaultdict

class ValidRule(BaseRule):
  NAME = "VALID"
  _RULES = BaseRule.RulesType()

  _QUALS2COPY="ID SRC TYPE SEQID START END STRAND LOCATION PHASE".split()

  @classmethod
  def prepare_context(cls, context):
    rules_data = context.get("_RULESDATA")
    if not rules_data:
      context.update({"_RULESDATA": defaultdict(dict)}, force_clean = True)
      rules_data = context.get("_RULESDATA")
    rules_data = rules_data["_ALL"]
    rules_data["USEDQUALS"] = None

  @classmethod
  def process(cls, context, re_context = None):
    # add stats
    context.global_context.add(cls.NAME, context)

  @classmethod
  def run_postponed(cls, context):
    for ctx in context.prev:
      if ctx.get("_ISLEAF"):
        context.used_leaves(ctx)

      used_quals = ctx.get("_RULESDATA")["_ALL"].get("USEDQUALS")
      if used_quals is None:
        ctx.get("_RULESDATA")["_ALL"]["USEDQUALS"] = {}
        used_quals = ctx.get("_RULESDATA")["_ALL"].get("USEDQUALS")

      for name in cls._QUALS2COPY:
        name = "_"+name
        value = ctx.get(name)
        if value:
          used_quals.update({name.lower():(name, value)})

      gff_quals = ctx.get("_QUALS")
      if gff_quals:
        for name, value in gff_quals.items():
          if name.lower() not in used_quals:
            if value:
              used_quals.update({name.lower():(name, value)})

    return


class ValidIfRule(ValidRule):
  NAME = "VALID_IF"
  _RULES = BaseRule.RulesType()
  # store gene.id/mrna.id/_feature at global context for checking

