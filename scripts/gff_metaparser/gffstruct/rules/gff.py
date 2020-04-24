from .base import BaseRule

import sys

class GffRule(BaseRule):
  # IDs can be substituted, but can be  ommited only for leaves
  # parents can't be ommited from GFF, paresnt, phase and strand  will be preserved in GFF
  # though parents should be excluded from qulifiers, they'll be added automatically based on parent feature ID
  NAME = "GFF"
  _RULES = BaseRule.RulesType()
  _FORCE_SUB = False

  @classmethod
  def prepare_context(cls, context):
    context.update(
      {
        "_%s_USEDQUALS" % cls.NAME : None,
        "_%s_QUALSCOPYALL" % cls.NAME : None,
      },
      force_clean = True,
    )

  def prepare_actions(self):
    self._target_quals = None
    raw = [ x.strip() for x in " ".join(self._actions_raw).replace(",", " ").split() if x.strip() ]
    raw = list(frozenset(raw))
    if raw:
      self._target_quals = raw

  def process(self, context, re_context = None):
    _uqname = "_%s_USEDQUALS" % self.NAME
    used_quals = context.get(_uqname)
    if used_quals is None:
      context.update({ _uqname:{} })
      used_quals = context.get(_uqname)

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
    _uqname = "_%s_USEDQUALS" % cls.NAME
    used_quals = context.get(_uqname)
    if used_quals is None:
      return
    is_leaf = context.get("_ISLEAF")
    if "phase" not in used_quals:
      phase = context.get("_PHASE")
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
    if context.get("_%s_QUALSCOPYALL" % cls.NAME):
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
      gff_uq = ctx.get("_GFF_USEDQUALS")
      gff_sub_uq = ctx.get("_GFF_SUB_USEDQUALS")
      #if gff_uq:
      #  print("GFF", gff_uq, file=sys.stderr)
      #if gff_sub_uq:
      #  print("GFF_SUB", gff_sub_uq, file=sys.stderr)
      if gff_uq:
        if gff_sub_uq:
          for k in gff_sub_uq:
            gff_uq[k] = gff_sub_uq[k]
          del ctx["_GFF_SUB_USEDQUALS"]
      elif gff_sub_uq:
        ctx["_GFF_USEDQUALS"] = gff_sub_uq
        del ctx["_GFF_SUB_USEDQUALS"]
      if ctx.get("_GFF_USEDQUALS") and ctx.get("_ISLEAF"):
        context.used_leaves(ctx)

