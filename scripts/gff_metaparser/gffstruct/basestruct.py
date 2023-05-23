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


# base structures

import sys

from collections import defaultdict

# local
from .rules import UnseenRule


# BASE STRUCTURES class
class BaseStructures:
    KNOWN_RULES = []

    def __init__(self, config, conf_patch=None, rule_options=None):
        self.rule_names = [r.NAME.lower().strip() for r in self.KNOWN_RULES]
        self.rule_factories = {r.NAME.lower().strip(): r for r in self.KNOWN_RULES}
        self.rule_options = frozenset(rule_options or [])

        # prepare factories from KNOWN_RULES
        self.config = config
        self.conf_patch = conf_patch
        self.load_conf()

        # propagate aliases to rules
        aliasCls = "alias" in self.rule_factories and self.rule_factories["alias"] or None

        # filling patterns to match to
        self.const_patterns = defaultdict(list)
        self.regex_patterns = list()

        for name in self.rule_factories:
            factory = self.rule_factories[name]
            factory.mature_regex(aliasCls)
            const_patterns = factory.const_patterns()
            there_were = frozenset(const_patterns.keys()) & frozenset(self.const_patterns.keys())
            if there_were:
                print("already seen pattern in config for %s" % (", ".join(there_were)), file=sys.stderr)
            for pat, rules in const_patterns.items():
                self.const_patterns[pat] += rules
            # regex_pattersn are checked dynamically for multiple hits
            self.regex_patterns += factory.regex_patterns()

    def parse_conf_str(self, raw):
        (nocmt, *_) = raw.partition("#")
        nocmt = nocmt.strip()
        if not nocmt:
            return None, None, None
        pattern, name, *actions = nocmt.split(maxsplit=2)
        return pattern, name, actions

    def load_conf(self):
        if not self.config:
            return
        patches = defaultdict(lambda: defaultdict(dict))
        if self.conf_patch:
            with open(self.conf_patch) as pf:
                for lineno, raw in enumerate(pf, start=1):
                    pattern, name, actions = self.parse_conf_str(raw)
                    if pattern is not None and name is not None:
                        patches[name][pattern]["actions"] = actions
                        patches[name][pattern]["lineno"] = "patch:%s" % lineno

        for lineno, raw in enumerate(self.config, start=1):
            pattern, name, actions = self.parse_conf_str(raw)
            if pattern is None:
                continue

            raw_name = name
            name_parts = name.split(":")
            name = name_parts[0]
            if len(name_parts) == 3:
                rule_opt = name_parts[1]
                not_rule_opt = False
                if rule_opt[0] == "!":
                    rule_opt = rule_opt[1:]
                    not_rule_opt = True
                if (rule_opt not in self.rule_options) ^ not_rule_opt:
                    name = name_parts[2]
                    print(
                        "using failback rule %s (from %s) for %s at %s" % (name, raw_name, pattern, lineno),
                        file=sys.stderr,
                    )
                else:
                    print(
                        "using main rule %s (from %s) for %s at %s" % (name, raw_name, pattern, lineno),
                        file=sys.stderr,
                    )
            elif len(name_parts) != 1:
                print(
                    "wrong number of name parts for rule %s at %s. skipping" % (raw_name, lineno),
                    file=sys.stderr,
                )
                continue

            if name in patches and pattern in patches[name]:
                actions = patches[name][pattern]["actions"]
                lineno_cfg = lineno
                lineno = patches[name][pattern]["lineno"]
                if actions == ["DISCARD"]:
                    print(
                        "discarding rule %s for %s at %s because of the patch DISCARD at %s"
                        % (name, pattern, lineno_cfg, lineno),
                        file=sys.stderr,
                    )
                    continue
                print(
                    "replacing rule %s for %s at %s because of the patch at %s"
                    % (name, pattern, lineno_cfg, lineno),
                    file=sys.stderr,
                )
            self.add_rule(name, pattern, actions, lineno)

    def add_rule(self, name_raw, pattern, actions, lineno):
        name = name_raw.lower().strip()
        if name not in self.rule_factories:
            print(
                "can't load unknown rule %s (raw: %s) at line %s" % (name, name_raw, lineno), file=sys.stderr
            )
            return
        factory = self.rule_factories[name]
        rule = factory(pattern, actions, lineno)

    def process(self, context, ignore_unseen=False):
        tag_raw = context.tag()
        tag = tag_raw.lower().strip()

        matched_rules = []
        if tag in self.const_patterns:
            matched_rules += list(map(lambda x: MatchedRuleCtx(x), self.const_patterns[tag]))
        # or stop if there's a constant match
        for it in self.regex_patterns:
            # print("matching regexp ", it.pattern, "tag", tag, file = sys.stderr )
            matching = it.re.fullmatch(tag_raw)
            if matching:
                matched_rules.append(MatchedRuleCtx(it.rule, matching))

        processed_rules = defaultdict(list)

        if matched_rules:
            for wrp in matched_rules:
                wrp.rule.process(context, re_context=wrp.re_context)
                processed_rules[tag].append(wrp.rule.NAME)
        elif not ignore_unseen:
            UnseenRule.process(context)
            processed_rules[tag].append(UnseenRule.NAME)

        return processed_rules

    def prepare_context(self, context):
        for name in self.rule_names:
            factory = self.rule_factories[name]
            factory.prepare_context(context)

    def prepare_postponed(self, context):
        for name in self.rule_names:
            factory = self.rule_factories[name]
            factory.prepare_postponed(context)

    def run_postponed(self, context):
        for name in self.rule_names:
            factory = self.rule_factories[name]
            # rule should decide itself if to run or not
            factory.run_postponed(context)


class MatchedRuleCtx:
    def __init__(self, rule, re_context=None):
        self.rule = rule
        self.re_context = re_context
