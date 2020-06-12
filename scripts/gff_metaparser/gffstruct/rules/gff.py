from .base import BaseRule

import sys
from collections import defaultdict

class GffRule(BaseRule):
  # IDs can be substituted, but can be  ommited only for leaves
  # parents can't be ommited from GFF, paresnt, phase and strand  will be preserved in GFF
  # though parents should be excluded from qulifiers, they'll be added automatically based on parent feature ID
  NAME = "GFF"
  _RULES = BaseRule.RulesType()
  _FORCE_SUB = False

  @classmethod
  def prepare_context(cls, context):
    rules_data = context.get("_RULESDATA")
    if not rules_data:
      context.update({"_RULESDATA": defaultdict(dict)}, force_clean = True)
      rules_data = context.get("_RULESDATA")
    rules_data = rules_data[cls.NAME]
    rules_data["USEDQUALS"] = None
    rules_data["QUALSCOPYALL"] = None


  def prepare_actions(self):
    self._target_quals = None
    raw = [ x.strip() for x in " ".join(self._actions_raw).replace(",", " ").split() if x.strip() ]
    raw = list(frozenset(raw))
    if raw:
      self._target_quals = raw

  def process(self, context, re_context = None):
    used_quals = context.get("_RULESDATA")[self.NAME].get("USEDQUALS")
    if used_quals is None:
      context.get("_RULESDATA")[self.NAME]["USEDQUALS"] = {}
      used_quals = context.get("_RULESDATA")[self.NAME].get("USEDQUALS")

    qname = context.get("_QNAME")
    if qname is None:
      context.update(_QUALSCOPYALL = True)
      return

    value = context.get("_LEAFVALUE")

    if not self._target_quals:
      used_quals.update({qname.lower():(qname, value)})
      return

    for new_name in self._target_quals:
      if not self._FORCE_SUB and new_name.lower() in used_quals:
        continue
      used_quals.update({new_name.lower():(new_name, value)})

  @classmethod
  def prepare_postponed(cls, context):
    # get seen quals and construct new qual
    used_quals = context.get("_RULESDATA")[cls.NAME].get("USEDQUALS")
    if used_quals is None:
      return
    is_leaf = context.get("_ISLEAF")
    if "phase" not in used_quals:
      phase = context.get("_PHASE")
      if phase is not None:
        used_quals.update({"phase":("phase", phase)})
    if not is_leaf:
      if "id" not in used_quals:
        used_quals.update({"ID":("ID", context.get("_ID"))})
    if context.get("_RULESDATA")[cls.NAME].get("QUALSCOPYALL"):
      # TODO: copy, everything not used
      # do not copy parent
      pass
    # print("prepare possponed: ", cls.NAME, str(used_quals), context.get("_ID"), context.get("_TYPE"), file=sys.stderr)

  @classmethod
  def run_postponed(cls, context):
    return


class GffSubRule(GffRule):
  NAME = "GFF_SUB"
  _RULES = BaseRule.RulesType()
  _FORCE_SUB = True

  @classmethod
  def run_postponed(cls, context):
    for ctx in context.prev:
      gff_uq = ctx["_RULESDATA"][GffRule.NAME].get("USEDQUALS")
      gff_sub_uq = ctx["_RULESDATA"][GffSubRule.NAME].get("USEDQUALS")

      # !!! no other rules can be used with the GFF or GFF_SUB rules !!!
      if gff_uq:
        ctx["_RULESDATA"]["_ALL"]["USEDQUALS"] = gff_uq
        del ctx["_RULESDATA"][GffRule.NAME]["USEDQUALS"]
        if gff_sub_uq:
          for k in gff_sub_uq:
            gff_uq[k] = gff_sub_uq[k]
          del ctx["_RULESDATA"][GffSubRule.NAME]["USEDQUALS"]
      elif gff_sub_uq:
        ctx["_RULESDATA"]["_ALL"]["USEDQUALS"] = gff_sub_uq
        del ctx["_RULESDATA"][GffSubRule.NAME]["USEDQUALS"]

      # update used_leaves
      if ctx["_RULESDATA"]["_ALL"].get("USEDQUALS") and ctx.get("_ISLEAF"):
        context.used_leaves(ctx)

