# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from .base import BaseRule

import re
import sys
from collections import defaultdict


class GffRule(BaseRule):
    # IDs can be substituted, but can be  ommited only for leaves
    # parents can't be ommited from GFF, paresnt, phase and strand  will be preserved in GFF
    # though parents should be excluded from qulifiers, they'll be added automatically based on parent feature ID
    NAME = "GFF"
    _RULES = BaseRule.RulesType()
    _FORCE_SUB = False
    _SPECIAL_KEYS = frozenset(["_FROM_STASH"])

    @classmethod
    def prepare_context(cls, context):
        rules_data = context.get("_RULESDATA")
        if not rules_data:
            context.update({"_RULESDATA": defaultdict(dict)}, force_clean=True)
            rules_data = context.get("_RULESDATA")
        rules_data = rules_data[cls.NAME]
        rules_data["USEDQUALS"] = None
        rules_data["QUALSCOPYALL"] = None

    def prepare_actions(self):
        self._target_quals = None

        raw = [x.strip() for x in "\t".join(self._actions_raw).split("\t") if x]
        if len(raw) < 1:
            return

        targets, *jsonstr = raw
        targets = [x.strip() for x in targets.replace(",", " ").split() if x.strip()]
        targets = list(frozenset(targets))
        if targets:
            self._target_quals = targets

        if not jsonstr:
            return

        addon, tech = self.parse_json(*jsonstr)
        self._actions = {
            "addon": addon,
            "from_stash": None,
        }
        self.add_actions_from_stash(tech)

    def add_actions_from_stash(self, tech):
        # i.e cds/parent GFF parent {"_FROM_STASH":[{"ID": "mrna:_PARENT/_STASH/cds_id"}]}
        if not tech or "_FROM_STASH" not in tech:
            return
        unstash = tech["_FROM_STASH"]
        out = []
        for qual, gctx_path in unstash.items():
            # special case: _IGNORE rule
            if qual == "_IGNORE":
                pat = list(filter(None, list(gctx_path)))
                if pat:
                    pat = "|".join(pat)
                    out.append({"qual": "_IGNORE", "ignore": re.compile(pat, flags=re.I)})
                    # print(f"adding ignore rule: '{pat}'", file=sys.stderr)
                continue
            # general path procesing
            pre_obj, *_path = gctx_path.split("/", 1)
            # if not path, use pre_obj as source ? todo
            obj_tag, *id_key = pre_obj.split(":")
            get_from_parent = False
            if not id_key:
                id_key = "_ID"
            elif id_key[0] == "_PARENT":
                id_key = pre_obj
                get_from_parent = True
            else:
                id_key = id_key[0]
            out.append(
                {
                    "qual": qual,
                    "obj_tag": obj_tag,
                    "from_parent": get_from_parent,
                    "id": id_key,
                    "path": _path and _path[0] or None,
                }
            )
        if out:
            self._actions["from_stash"] = out

    def process(self, context, re_context=None):
        used_quals = context.get("_RULESDATA")[self.NAME].get("USEDQUALS")
        if used_quals is None:
            context.get("_RULESDATA")[self.NAME]["USEDQUALS"] = {}
            used_quals = context.get("_RULESDATA")[self.NAME].get("USEDQUALS")

        qname = context.get("_QNAME")
        if qname is None:
            context.update(_QUALSCOPYALL=True)
            return

        value = context.get("_LEAFVALUE")

        if not self._target_quals and not self._actions:
            used_quals.update({qname.lower(): (qname, value)})
            return

        ignore_re = None
        if self._actions:
            from_stash = self._actions.get("from_stash", [])
            # get ignore pattern first (only first)
            ignore_list = list(filter(lambda x: x.get("qual") == "_IGNORE" and x.get("ignore"), from_stash))
            if ignore_list:
                ignore_re = ignore_list[0].get("ignore")
            # process qual list
            for it in from_stash:
                name = it["qual"]
                if name == "_IGNORE":
                    continue
                value = self.render_stashed(it, context)
                if not self._FORCE_SUB and name.lower() in used_quals:
                    continue
                if ignore_re:
                    value = list(
                        filter(lambda x: not isinstance(x, str) or not ignore_re.search(x), list(value))
                    )
                if not value:
                    continue
                if value:
                    used_quals.update({name.lower(): (name, value)})

        for new_name in self._target_quals:
            if not self._FORCE_SUB and new_name.lower() in used_quals:
                continue
            # ignore pattern check
            if ignore_re:
                # print(f"match {value} against {ignore_re}", file=sys.stderr)
                value = list(filter(lambda x: not isinstance(x, str) or not ignore_re.search(x), list(value)))
                # print(f"matched {value}", file=sys.stderr)
            if not value:
                continue
            used_quals.update({new_name.lower(): (new_name, value)})

    def render_stashed(self, stashed, context=None):
        if not stashed or not context:
            return None
        _id = "%s%s" % (self._CTX_PFX, stashed.get("id"))
        _id = context.get(_id)
        if not _id:
            return None
        # try to get data from global context
        val = context.global_context.get(stashed["obj_tag"], _id, stashed["path"])
        if val is not None:
            return val
        out = stashed.copy()
        out["id"] = _id
        return {"_FROM_STASH": out}

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
                used_quals.update({"phase": ("phase", phase)})
        if not is_leaf:
            if "id" not in used_quals:
                used_quals.update({"ID": ("ID", context.get("_ID"))})
        if context.get("_RULESDATA")[cls.NAME].get("QUALSCOPYALL"):
            # TODO: copy, everything not used
            # do not copy parent
            pass
        # print("prepare postponed: ", cls.NAME, str(used_quals), context.get("_ID"), context.get("_TYPE"), file=sys.stderr)

    @classmethod
    def run_postponed(clsf, context, name_override=None):
        return


class GffSubRule(GffRule):
    NAME = "GFF_SUB"
    _RULES = BaseRule.RulesType()
    _FORCE_SUB = True

    @classmethod
    def run_postponed(cls, context, name_override=None):
        # id_tag = "%s%s" % (self._CTX_PFX, id_tag)
        # obj_id = context.get(id_tag)
        # print(context.data.keys(), file = sys.stderr)
        # use from stash / store before in process
        # ? update from global ctx in gff write, because only there there's filled global contex. though, prepare in process

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
