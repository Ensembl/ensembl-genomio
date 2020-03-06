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
    # use self.rule_factories['alias'].alias2regexp()

    # filling patterns to match to
    self.const_patterns = defaultdict()
    self.regex_patterns = []
    for name in self.rule_factories:
      factory = self.rule_factories[name]
      const_patterns = factory.const_patterns()
      rewriting = frozenset(const_patterns) & frozenset(self.const_patterns)
      if rewriting:
         print("const_patterns clash for %s" % (" ,".join(rewriting)) ,file=sys.stderr)
      else:
        self.const_patterns.update(const_patterns)
      # check regex? or dynamically when multiple hits
      regex_patterns = factory.regex_patterns
      if factory.regex_patterns:
        self.regex_patterns.append(regex_patterns)

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
    factory = self.rule_factories[name]
    rule = factory(pattern, actions, lineno)

    rule.dump()


  def process(self, struct):
    tag = struct.tag
    rule, re_context = None, None
    if tag in self.const_patterns:
      rule = self.const_patterns[tag]
    else:
      # do not process if there was a static match, match first
      for it in self.regex_patterns:
        re_context = it.regex_patterns(tag)
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
