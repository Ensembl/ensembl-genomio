from .base import BaseRule

import sys

class GffRule(BaseRule):
  # IDs can be substituted, but can be  ommited only for leaves
  # parents can't be ommited from GFF, paresnt, phase and strand  will be preserved in GFF
  # though parents should be excluded from qulifiers, they'll be added automatically based on parent feature ID
  NAME = "GFF"
  _RULES = BaseRule.RulesType()

  @classmethod
  def prepare_context(cls, context):
    context.update(
      force_clean = True,
      _USEDQUALS = None,
      _QUALSCOPYALL = None,
    )

  def prepare_actions(self):
     pass

  def process(self, context, re_context = None):
   used_quals = context.get("_USEDQUALS")
   if used_quals is None:
     context.update(_USEDQUALS = {})
     used_quals = context.get("_USEDQUALS")
   qname = context.get("_QNAME")
   if qname is None:
     context.update(_QUALSCOPYALL = True)
     return
   #new_name = self._new_name
   # TODO: use new name, what if changing feature?
   used_quals.update({qname.lower():(qname, context.get("_LEAFVALUE"))})

  @classmethod
  def prepare_postponed(cls, context):
    # get seen quals and construct new qual
    used_quals = context.get("_USEDQUALS")
    if used_quals is None:
      return
    is_leaf = context.get("_ISLEAF")
    if "phase" not in used_quals:
      phase = context.get("_PHASE")
      if phase is not None:
        used_quals.update({"phase":phase})
    if not is_leaf:
      if "id" not in used_quals:
        used_quals.update({"ID":context.get("_ID")})
    if context.get("_QUALSCOPYALL"):
      # TODO: copy, everything not used
      # do not copy parent
      pass
    print("prepare possponed: ", str(used_quals), context.get("_ID"), context.get("_TYPE"), file=sys.stderr)


class GffSubRule(GffRule):
  NAME = "GFF_SUB"
  _RULES = BaseRule.RulesType()

  @classmethod
  def prepare_postponed(cls, context):
    pass

