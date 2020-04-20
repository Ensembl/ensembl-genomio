# base structures

import sys

from collections import defaultdict

# local
from .rules import UnseenRule

# BASE STRUCTURES class
class BaseStructures:
  KNOWN_RULES = [
  ]

  def __init__(self, config, conf_patch = None):
    self.rule_names = [ r.NAME.lower().strip() for r in self.KNOWN_RULES ]
    self.rule_factories = { r.NAME.lower().strip() : r for r in self.KNOWN_RULES }

    # prepare factories from KNOWN_RULES
    self.config = config
    self.conf_patch = conf_patch
    self.load_conf()

    # propagate aliases to rules
    aliasCls = 'alias' in self.rule_factories and self.rule_factories['alias'] or None

    # filling patterns to match to
    self.const_patterns = defaultdict(list)
    self.regex_patterns = list()

    for name in self.rule_factories:
      factory = self.rule_factories[name]
      factory.mature_regex(aliasCls)
      const_patterns = factory.const_patterns()
      there_were = frozenset(const_patterns.keys()) & frozenset(self.const_patterns.keys())
      if there_were:
         print("already seen pattern in config for %s" % (", ".join(there_were)) ,file=sys.stderr)
      for pat, rules in const_patterns.items():
        self.const_patterns[pat] += rules
      # regex_pattersn are checked dynamically for multiple hits
      self.regex_patterns += factory.regex_patterns()

  def parse_conf_str(self, raw):
    (nocmt, *_) = raw.partition('#')
    nocmt = nocmt.strip()
    if not nocmt:
      return None, None, None
    pattern, name, *actions = nocmt.split(maxsplit=2)
    return pattern, name, actions

  def load_conf(self):
    if not self.config:
      return
    patches = defaultdict(lambda:defaultdict(dict))
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
      if name in patches and pattern in patches[name]:
        actions = patches[name][pattern]["actions"]
        lineno = patches[name][pattern]["lineno"]
      self.add_rule(name, pattern, actions, lineno)

  def add_rule(self, name_raw, pattern, actions, lineno):
    name = name_raw.lower().strip()
    if name not in self.rule_factories:
      print("can't load unknown rule %s (raw: %s) at line %s" % (name, name_raw, lineno), file = sys.stderr)
      return
    factory = self.rule_factories[name]
    rule = factory(pattern, actions, lineno)

  def process(self, context, ignore_unseen = False):
    tag_raw = context.tag()
    tag = tag_raw.lower().strip()

    matched_rules = []
    if tag in self.const_patterns:
      matched_rules += list(map(lambda x: MatchedRuleCtx(x), self.const_patterns[tag]))
    # or stop if there's a constant match
    for it in self.regex_patterns:
      matching = it.re.fullmatch(tag)
      if matching:
        matched_rules.append(MatchedRuleCtx(it.rule, matching))

    processed_rules = defaultdict(list)

    if matched_rules:
      for wrp in matched_rules:
        wrp.rule.process(context, re_context = wrp.re_context)
        processed_rules[tag].append(wrp.rule.NAME)
    elif not ignore_unseen:
      UnseenRule.process(context, noconfig = (self.config == None))
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
  def __init__(self, rule, re_context = None):
    self.rule = rule
    self.re_context = re_context

