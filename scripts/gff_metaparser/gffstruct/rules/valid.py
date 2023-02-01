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


import sys

from collections import defaultdict
from .base import BaseRule


class ValidRule(BaseRule):
    NAME = "VALID"
    _RULES = BaseRule.RulesType()

    _QUALS2COPY = "ID SRC TYPE SEQID START END STRAND LOCATION PHASE".split()

    @classmethod
    def prepare_context(cls, context):
        rules_data = context.get("_RULESDATA")
        if not rules_data:
            context.update({"_RULESDATA": defaultdict(dict)}, force_clean=True)
            rules_data = context.get("_RULESDATA")
        rules_data = rules_data[cls.NAME]
        rules_data["USEDQUALS"] = None

    @classmethod
    def process(cls, context, re_context=None, name_override=None):
        _name = name_override or cls.NAME
        # add stats
        context.global_context.add(_name, context)
        context.run_to_root(updater=lambda x: cls.mark_as_usefull(x, name_override=_name))
        context.update(force_clean=True, _RECTX=re_context)

    @classmethod
    def mark_as_usefull(cls, ctx, name_override=None):
        _name = name_override or cls.NAME
        used_quals = ctx.get("_RULESDATA")[_name].get("USEDQUALS")
        if used_quals is None:
            ctx.get("_RULESDATA")[_name]["USEDQUALS"] = {}

    @classmethod
    def run_postponed(cls, context, name_override=None):
        name_to_check = name_override or cls.NAME
        for ctx in context.prev:
            check_quals = ctx.get("_RULESDATA")[name_to_check].get("USEDQUALS")
            if check_quals is None:
                continue

            if ctx.get("_ISLEAF"):
                context.used_leaves(ctx)

            used_quals = ctx.get("_RULESDATA")["_ALL"].get("USEDQUALS")
            if used_quals is None:
                ctx.get("_RULESDATA")["_ALL"]["USEDQUALS"] = {}
                used_quals = ctx.get("_RULESDATA")["_ALL"].get("USEDQUALS")

            for name in cls._QUALS2COPY:
                name = "_" + name
                value = ctx.get(name)
                if value:
                    used_quals.update({name.lower(): (name, value)})

            gff_quals = ctx.get("_QUALS")
            if gff_quals:
                for name, value in gff_quals.items():
                    if name.lower() not in used_quals:
                        if value:
                            used_quals.update({name.lower(): (name, value)})

        return


class ValidIfRule(ValidRule):
    NAME = "VALID_IF"
    _RULES = BaseRule.RulesType()
    # store gene.id/mrna.id/_feature at global context for checking
