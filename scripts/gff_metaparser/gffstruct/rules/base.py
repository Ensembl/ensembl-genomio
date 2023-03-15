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


# rules classes

import json
import re
import sys

from collections import defaultdict
from copy import deepcopy, copy


def _RulesType():
    return {
        "const_match": defaultdict(list),
        "_regex_match_pre": defaultdict(list),
        "regex_match": list(),
    }


class BaseRule:
    NAME = ""
    _RULES = _RulesType()
    _SPECIAL_KEYS = []
    _CTX_PFX = "_TECH_4KIDS_"

    @classmethod
    def RulesType(cls):
        return _RulesType()

    @classmethod
    def AliasSymbol(cls):
        return r"@"

    def __init__(self, pattern, actions, lineno):
        self._pattern = pattern
        self._actions_raw = actions
        self._actions = None
        self._lineno = lineno
        self.update_rules()
        self.prepare_actions()

    def update_rules(self):
        pat = self._pattern.strip().lower()
        rules = self._RULES["const_match"]
        if BaseRule.AliasSymbol() in pat:
            rules = self._RULES["_regex_match_pre"]

        if pat in rules:
            print(
                "already seen pattern %s at lines %s for %s at %d"
                % (
                    self._pattern,
                    ",".join(map(lambda r: str(r._lineno), rules[pat])),
                    self.NAME,
                    self._lineno,
                ),
                file=sys.stderr,
            )
        rules[pat].append(self)

    @classmethod
    def mature_regex(cls, aliases_cls):
        if not aliases_cls:
            return
        for pat, rules in cls._RULES["_regex_match_pre"].items():
            for rule in rules:
                raw_pat = rule._pattern.strip()

                re_pat = aliases_cls.mature_regex(raw_pat)
                # print("maturing pat ", pat, " raw_pat ", raw_pat, " re_pat ", re_pat, file=sys.stderr)

                if re_pat == None:
                    continue
                if re_pat == raw_pat:
                    cls._RULES["const_match"][pat].append(rule)
                    print(
                        r"cannot mature pattern %s (raw: %s) for %s (line %d)"
                        % (
                            pat,
                            raw_pat,
                            cls.NAME,
                            rule._lineno,
                        ),
                        file=sys.stderr,
                    )
                else:
                    cls._RULES["regex_match"].append(RERuleWrapper(rule, re_pat))

    @classmethod
    def const_patterns(cls):
        return copy(cls._RULES["const_match"])

    @classmethod
    def regex_patterns(cls):
        return copy(cls._RULES["regex_match"])

    def prepare_actions(self):
        pass

    def process(self, context, re_context=None):
        print(
            "processing %s for %s with match groups %s"
            % (context.tag(), self.NAME, str(re_context and re_context.groupdict() or None))
        )
        pass

    @classmethod
    def prepare_context(cls, context):
        pass

    @classmethod
    def prepare_postponed(cls, context):
        pass

    @classmethod
    def run_postponed(clsf, context, name_override=None):
        pass

    def parse_json(self, *s):
        data, tech = None, None
        if not s:
            return data, tech

        try:
            raw = json.loads(s[0])
            data = {k: v for k, v in raw.items() if k not in self._SPECIAL_KEYS}
            tech = {k: v for k, v in raw.items() if k in self._SPECIAL_KEYS}
        except:
            raise Exception("wrong JSON part for rule %s at line %s: %s" % (self.NAME, self._lineno, str(s)))

        return data, tech


# ideally, add wrapping for static, as well
class RuleWrapper:
    def __init__(self, rule, pattern=None, const_match=True):
        self.rule = rule
        self.pattern = pattern
        self.const_match = const_match


class RERuleWrapper(RuleWrapper):
    def __init__(self, rule, pattern=None, const_match=False):
        self.rule = rule
        self.pattern = pattern
        self.const_match = const_match
        self.re = None
        if pattern:
            self.re = re.compile(pattern, flags=re.I)
