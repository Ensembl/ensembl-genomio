from .base import BaseRule

import json
import re
import sys

class JsonRule(BaseRule):
  NAME = "JSON"
  _RULES = BaseRule.RulesType()
  _SPECIAL_KEYS =  frozenset(["_IGNORE","_MAP", "_S"])
  _OUTPUT_FORCE_SUB = False

  def parse_json(self, *s):
    data, tech = None, None
    if not s:
      return data, tech

    try:
      raw = json.loads(s[0])
      data = { k: v for k,v in raw.items() if k not in self._SPECIAL_KEYS }
      tech = { k: v for k,v in raw.items() if k in self._SPECIAL_KEYS }
    except:
      raise Exception("wrong JSON part for rule %s at line %s: %s" % (self.NAME, self._lineno, str(s)))

    return data, tech

  def prepare_actions(self):
    raw = [ x.strip() for x in "\t".join(self._actions_raw).split("\t") if x ]
    if len(raw) < 1:
      print("no actions for rule %s at line %s" % (self.NAME, self._lineno), file = sys.stderr)
      return
    tag_path, *jsonstr = raw
    addon, tech = self.parse_json(*jsonstr)

    if tag_path.count("@")  != tag_path.count("/@"):
      raise Exception("wrong @list notation used in path string for rule %s at line %s" % (self.NAME, self._lineno))
    if tag_path.count("@") > 1:
      raise Exception("to many @list in path string for rule %s at line %s" % (self.NAME, self._lineno))

    out_tag_raw, *path = tag_path.split("/", 1)
    if not path:
      path = ""
    else:
      path = path[0]
    # id part of the out_tag 
    store_to_parent = False
    out_tag, *id_key = out_tag_raw.split(":")
    if not id_key:
      id_key = "_ID"
    elif id_key[0] == "_PARENT":
      id_key = out_tag_raw
      store_to_parent = True
    else:
      id_key = id_key[0]

    self._actions = {
      "out_tag" : out_tag,
      "store_to_parent" : store_to_parent,
      "id_key" : id_key,
      "path" : path,
      "addon" : addon, # interpolate addon, run sub, run map afterwards
      "map" : tech and tech.get("_MAP", None) or None
    }
    self.add_actions_sub(tech)
    self.add_actions_ignore(tech)

  def add_actions_ignore(self, tech):
    if not tech or "_IGNORE" not in tech:
      return
    ign = tech["_IGNORE"]
    if not ign:
      return
    if type(ign) == list:
      ign = "|".join(ign)
    self._actions["ignore"] = re.compile(ign, flags = re.I)

  def add_actions_sub(self, tech):
    if not tech or "_SUB" not in tech:
      return
    sub = tech["_SUB"]
    if not sub:
      return
    if type(sub) != list:
      sub = [sub]
    self._actions["sub"] = []
    for from_to in sub:
      from_, to_ = from_to.split(";", 1)
      from_ = re.compile(from_, flags = re.I)
      self._actions["sub"].append((from_, to_))

  def interpolate_values(self, d):
    pass

  def process(self, context, re_context = None):
    if not self._actions:
      print("no actions for rule %s at line %s applyed to %s" % (self.NAME, self._lineno, context.tag()), file = sys.stderr)
    return
    # get addon
    # add context._leafvalue to the addon if not None
    # interpolate JSON values
    # check ignore, check map
    # context.get("_JSON_OUTTAG_IDS", {})
    # context.add("_JSON_OUTTAG_IDS" = )
    # get object id (from context or from parent context)
    # context.global_context(obj_tag, obj_id, json, **global_context_add_flags)
    # context.global_context(obj_tag, obj_id, path, json, substitute = self._OUTPUT_FORCE_SUB)


class JsonSubRule(JsonRule):
  NAME = "JSON_SUB"
  _RULES = JsonRule.RulesType()
  _SPECIAL_KEYS = JsonRule._SPECIAL_KEYS # frozenset(["_IGNORE","_MAP", "_S"])
  _OUTPUT_FORCE_SUB = True


