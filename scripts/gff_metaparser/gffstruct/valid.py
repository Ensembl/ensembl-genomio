# valid structures

import sys

from collections import defaultdict

from .rules import *

# VALID STRUCTURES class

class ValidStructures:
  KNOWN_RULES = [
    AliasRule,
    IgnoreRule,
    ValidRule,
    ValidIfRule,
    FixRule,
    ForceFixRule,
    SubRule,
    ForceSubRule,
    SetRule,
  ]

  def __init__(self, config):
    self.rule_factories = { r.NAME.lower().strip() : r for r in self.KNOWN_RULES }

    # prepare factories from KNOWN_RULES
    self.config = config
    self.load_conf()

    # propagate aliases to rules

    # filling patterns to match to
    self.const_match = defaultdict()
    self.regex_match = []
    for name in self.rule_factories:
      rule = self.rule_factories[name]
      self.const_match.update([ (pat, rule) for pat in rule.const_patterns() ])
      if rule.regex_patterns:
        self.regex_match.append(rule)

  def load_conf(self):
    if not self.config:
      return
    for lineno, raw in enumerate(self.config, start=1):
      (nocmt, *_) = raw.partition('#')
      nocmt = nocmt.strip()
      if not nocmt:
        continue
      pattern, name, *actions = nocmt.split(maxsplit=2)
      self.add_rule(name, pattern, actions, lineno)

  def add_rule(self, name_raw, pattern, actions, lineno):
    name = name_raw.lower().strip()
    if name not in self.rule_factories:
       print("can't load unknown rule %s at line %s" % (name, lineno), file = sys.stderr)
       return
    fabric = self.rule_factories[name]
    rule = fabric(pattern, actions, lineno)

    rule.dump()


  def process(self, struct):
    tag = struct.tag
    rule, re_context = None, None
    if tag in self.const_match:
      rule = self.const_match[tag]
    else:
      # do not process if there was a static match, match first
      for it in self.regex_match:
        re_context = it.regex_match(tag)
        if re_context:
          rule = it
          break
    if rule:
      rule.process(struct, re_context = re_context)
    else:
      # log not seen stuct / struct tag
      UnseenRule.process(struct, noconfig = (self.config == None))
      pass
    # return ???
