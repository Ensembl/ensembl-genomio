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


#
import sys
import re

from .base import BaseRule

from collections import defaultdict


class AliasRule(BaseRule):
    NAME = "ALIAS"
    _RULES = BaseRule.RulesType()

    aliases = defaultdict(str)
    capture = dict()
    regexp = dict()
    lineno = defaultdict(list)

    @classmethod
    def update_aliases(cls, pattern, actions):
        pattern = pattern.strip()

        pre_pat = ",".join(actions).split(",")
        pre_pat = map(lambda x: x.strip(), pre_pat)
        pre_pat = filter(None, pre_pat)
        pre_pat = map(
            lambda x: r"\\b%s\\b" % x, filter(None, pre_pat)
        )  # N.B. double quoting, because used in re.sub
        pre_pat = list(pre_pat)

        if pattern in cls.aliases and pre_pat:
            cls.aliases[pattern] += "|"
        cls.aliases[pattern] += "|".join(pre_pat)

        capture = r"(?P<%s>%s)" % (pattern.replace(BaseRule.AliasSymbol(), ""), cls.aliases[pattern])
        cls.capture[pattern] = capture
        rx = r"%s\b" % pattern  # should start with '@' otherwise add initial \b
        # print("rx %s capture %s pattern %s " %(rx, capture, pattern), file=sys.stderr)
        cls.regexp[pattern] = re.compile(rx)

    def __init__(self, pattern, actions, lineno=None):
        if not pattern or not pattern.startswith("@"):
            print(
                r"ignoring: rule %s for %s (line %d): alias should start with '@'"
                % (self.NAME, pattern, lineno),
                file=sys.stderr,
            )
            return
        AliasRule.update_aliases(pattern, actions)
        AliasRule.update_lineno(pattern, lineno)

    @classmethod
    def update_lineno(cls, pattern, lineno):
        cls.lineno[pattern].append(lineno)

    @classmethod
    def mature_regex(cls, pattern):
        if pattern == AliasRule:
            return None
        # try case insensitive
        for alias, capture in cls.capture.items():
            rx = cls.regexp[alias]
            pattern = rx.sub(capture, pattern)
        return pattern
