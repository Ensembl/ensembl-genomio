from .base import BaseRule

import json
import re
import sys

class JsonRule(BaseRule):
  NAME = "JSON"
  _RULES = BaseRule.RulesType()
  _SPECIAL_KEYS =  frozenset(["_IGNORE","_MAP", "_S"])
  _OUTPUT_FORCE_SUB = False
  _CTX_PFX="_TECH_JSON_"

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


    key = ""
    out_tag_raw, *path = tag_path.split("/")
    if not path:
      path = ""
    else:
      *path, key = path
      path = "/".join(path)

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
      "key" : key,
      "addon" : addon, # interpolate addon, run sub, run map afterwards
      "map" : tech and tech.get("_MAP", None) or None,
      "sub" : None,
      "ignore" : None,
    }
    self.add_actions_sub(tech)
    self.add_actions_ignore(tech)

  def interpolate(self, data, context):
    interpolated = False
    for k in data:
      v = data[k]
      if v is None or type(v) != str:
        continue
      v = context.get(v, default = v)

      asub = self._actions["sub"]
      amap = self._actions["map"]
      aignore = self._actions["ignore"]
      if asub:
        for f,t in asub:
          v = f.sub(t, v)
      if amap and v in amap:
        v = amap[v]
      if aignore and aignore.search(v) is not None:
        v = None
      if v != data[k]:
        interpolated = True
      data[k] = v
    if interpolated:
      print(data, file = sys.stderr)
    return interpolated

  def process(self, context, re_context = None):
    if not self._actions:
      print("no actions for rule %s at line %s applyed to %s" % (self.NAME, self._lineno, context.tag()), file = sys.stderr)
      return

    _a = self._actions

    obj_tag = _a["out_tag"]
    id_tag = _a["id_key"]
    tag4kids = "%s%s:_PARENT" % (self._CTX_PFX, obj_tag) # only valid if not using parent :(

    obj_id = None
    if _a["store_to_parent"]:
      id_tag = "%s%s" % (self._CTX_PFX, id_tag)
      obj_id = context.get(id_tag)
    else:
      obj_id = context.get(id_tag)
      context.update({tag4kids:obj_id}) # store for kids

    value = context.get("_LEAFVALUE")
    if type(value) == list and len(value) == 1:
      value = value[0]

    data = {}
    if _a["key"]:
      data = { _a["key"] : value }
    if _a["addon"]:
      data.update(_a["addon"])
    if self.interpolate(data, context):
      print("interpolating data for rule %s at line %s applyed to %s" % (self.NAME, self._lineno, context.tag()), file = sys.stderr)

    context.global_context.add(obj_tag, obj_id, _a["path"], data, force = self._OUTPUT_FORCE_SUB)


class JsonSubRule(JsonRule):
  NAME = "JSON_SUB"
  _RULES = JsonRule.RulesType()
  _OUTPUT_FORCE_SUB = True


