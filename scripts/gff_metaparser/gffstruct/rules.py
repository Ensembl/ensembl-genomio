# rules classes
from collections import defaultdict

class BaseRule:
  NAME = ""
  _const_patterns = {}
  _regex_patterns = []

  def __init__(self, pattern, actions, lineno):
    pass

  @classmethod
  def const_patterns(cls):
    return cls._const_patterns.keys()

  @classmethod
  def regex_patterns(cls):
    return list(map(lambda x: x[0], cls._regex_patterns))

  def const_match(self, tag):
    pass

  def regex_match(self, tag):
    pass

  def process(self, tag, re_context = None):
    pass

  @classmethod
  def process_unknown(cls, struct, noconfig=True):
    # log not seen stuct / struct tag
    pass

  def dump(self):
    pass


class AliasRule(BaseRule):
  NAME = "ALIAS"
  aliases = defaultdict(str)
  lineno = defaultdict(list)

  @classmethod
  def update_aliases(cls, pattern, actions):
    pattern = pattern.strip()

    pre_pat = ",".join(actions).split(",")
    pre_pat = map(lambda x: x.strip(), pre_pat)
    pre_pat = filter(None, pre_pat)
    pre_pat = map(lambda x: r'\b%s\b' % x, filter(None, pre_pat))

    cls.aliases[pattern] = "|".join( [cls.aliases[pattern]] + list(pre_pat) )


  def __init__(self, pattern, actions, lineno = None):
     if not pattern or not pattern.startswith("@"):
       print("ignoring: rule %s for %s (line %d): alias should start with '@'" % (
              self.NAME, pattern, lineno
            ), file=sys.stderr)
       return
     AliasRule.update_aliases(pattern, actions)
     AliasRule.update_lineno(pattern, lineno)

  @classmethod
  def update_lineno(cls, pattern, lineno):
    cls.lineno[pattern].append(lineno)

  def dump(self):
    print(self.aliases)
    print(self.lineno)


class IgnoreRule(BaseRule):
  NAME = "IGNORE"


