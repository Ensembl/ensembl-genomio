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


from .base import BaseRule

import re
import sys


class JsonRule(BaseRule):
    NAME = "JSON"
    _RULES = BaseRule.RulesType()
    _SPECIAL_KEYS = frozenset(["_IGNORE", "_MAP", "_SUB", "_SPLIT", "_NUMVAL", "_PREPEND", "_APPEND"])
    _OUTPUT_FORCE_SUB = False

    def add_actions_ignore(self, tech):
        if not tech or "_IGNORE" not in tech:
            return
        ign = tech["_IGNORE"]
        if not ign:
            return
        if type(ign) == list:
            ign = "|".join(ign)
        self._actions["ignore"] = re.compile(ign, flags=re.I)

    def add_actions_split(self, tech):
        if not tech or "_SPLIT" not in tech:
            return
        split = tech["_SPLIT"]
        if not split:
            return
        delim, *keys = split.split(";")
        self._actions["split"] = {"delim": delim, "keys": keys}

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
            from_ = re.compile(from_, flags=re.I)
            self._actions["sub"].append((from_, to_))

    @classmethod
    def a2n(cls, val, totype=None):
        if val is None or type(val) != str:
            return val
        if not totype:
            return val
        if totype != True and totype != "int":
            return val
        try:
            return int(val)
        except:
            pass
        try:
            outval = float(val)
            if totype == "int":
                return int(outval)
            return outval
        except:
            pass
        return val

    def add_actions_numval(self, tech):
        if not tech or "_NUMVAL" not in tech:
            return
        numval = tech["_NUMVAL"]
        if not numval:
            return
        self._actions["numval"] = lambda x: JsonRule.a2n(x, totype=numval)

    def add_actions_prepend(self, tech):
        if not tech or "_PREPEND" not in tech:
            return
        prefix = tech["_PREPEND"]
        if not prefix:
            return
        self._actions["prepend"] = {"prefix": prefix}

    def add_actions_append(self, tech):
        if not tech or "_APPEND" not in tech:
            return
        suffix = tech["_APPEND"]
        if not suffix:
            return
        self._actions["append"] = {"suffix": suffix}

    def prepare_actions(self):
        raw = [x.strip() for x in "\t".join(self._actions_raw).split("\t") if x]
        if len(raw) < 1:
            print("no actions for rule %s at line %s" % (self.NAME, self._lineno), file=sys.stderr)
            return
        tag_path, *jsonstr = raw
        addon, tech = self.parse_json(*jsonstr)

        if tag_path.count("@") != tag_path.count("/@"):
            raise Exception(
                "wrong @list notation is used in path string for rule %s at line %s"
                % (self.NAME, self._lineno)
            )
        if tag_path.count("@") > 1:
            raise Exception("to many @list in path string for rule %s at line %s" % (self.NAME, self._lineno))
        if tag_path.count("!") != tag_path.count("/!"):
            raise Exception(
                "wrong !scalar notation is used in path string for rule %s at line %s"
                % (self.NAME, self._lineno)
            )
        if tag_path.count("!") > 1:
            raise Exception(
                "to many !scalar in path string for rule %s at line %s" % (self.NAME, self._lineno)
            )
        if tag_path.count("!") == 1 and not tag_path.split("/")[-1].startswith("!"):
            raise Exception(
                "!scalar not as the last part of path string for rule %s at line %s"
                % (self.NAME, self._lineno)
            )

        key = ""
        out_tag_raw, *path = tag_path.split("/")
        if not path:
            path = ""
        else:
            *path, key = path
            path = "/".join(path)

        key_is_list = None
        if key.startswith("@"):  # force list
            key = key[1:]
            key_is_list = True
        if key.startswith("!"):  # force scalar
            key = key[1:]
            key_is_list = False

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
            "out_tag": out_tag,
            "store_to_parent": store_to_parent,
            "id_key": id_key,
            "path": path,
            "key": key,
            "key_is_list": key_is_list,
            "addon": addon,  # interpolate addon, run sub, run map afterwards
            "map": tech and tech.get("_MAP", None) or None,
            "sub": None,
            "ignore": None,
            "split": None,
            "numval": None,
            "append": None,
            "prepend": None,
        }
        self.add_actions_prepend(tech)
        self.add_actions_append(tech)
        self.add_actions_split(tech)
        self.add_actions_sub(tech)
        self.add_actions_ignore(tech)
        self.add_actions_numval(tech)

    def interpolate(self, data, context, do_split_and_earlier_actions=True):
        interpolated = False

        if data is None:
            return None, interpolated

        if type(data) == dict:
            out = dict()
            for k in data:
                v, _interpolated = self.interpolate(data[k], context, do_split_and_earlier_actions)
                interpolated = interpolated or _interpolated
                if v is not None:
                    out[k] = v
                else:
                    # drop the whole dict, if field is gone
                    if data[k] is not None:
                        return None, True
            if not out:
                out = None
            return out, interpolated

        if type(data) == list:
            out = []
            for v in data:
                v, _interpolated = self.interpolate(v, context, do_split_and_earlier_actions)
                interpolated = interpolated or _interpolated
                if v is not None:
                    out.append(v)
            if not out:
                out = None
            return out, interpolated

        # scalar section
        v = context.get(data, data)
        # changing from 0-based internal (start-1, end] to 1-based external for JSON [start, end]
        if data == "_START" and data != v:
            v += 1

        aprepend = self._actions["prepend"]
        aappend = self._actions["append"]
        asplit = self._actions["split"]
        asub = self._actions["sub"]
        amap = self._actions["map"]
        aignore = self._actions["ignore"]
        anumval = self._actions["numval"]

        if do_split_and_earlier_actions and type(v) == str and aprepend and aprepend.get("prefix"):
            v = str(aprepend["prefix"]) + v

        if do_split_and_earlier_actions and type(v) == str and aappend and aappend.get("suffix", ""):
            v = v + str(aappend[suffix])

        if asplit and do_split_and_earlier_actions:
            res = v.split(asplit["delim"])
            if asplit["keys"]:
                return self.interpolate(
                    dict(zip(asplit["keys"], res)), context, do_split_and_earlier_actions=False
                )
            else:
                return self.interpolate(res, context, do_split_and_earlier_actions=False)

        if asub and v and type(v) == str:
            for f, t in asub:
                if v:
                    v = f.sub(t, v)
        if amap:
            if v in amap:
                v = amap[v]
            elif amap.get("_IGNORE_REST"):
                v = None
        if aignore and v and aignore.search(v) is not None:
            v = None
        if anumval:
            v = anumval(v)

        return v, v != data

    def process(self, context, re_context=None):
        if not self._actions:
            print(
                "no actions for rule %s at line %s applyed to %s" % (self.NAME, self._lineno, context.tag()),
                file=sys.stderr,
            )
            return

        _a = self._actions

        obj_tag = _a["out_tag"]
        id_tag = _a["id_key"]
        tag4kids = "%s%s:_PARENT" % (self._CTX_PFX, obj_tag)  # only valid if not using parent :(

        obj_id = None
        if _a["store_to_parent"]:
            id_tag = "%s%s" % (self._CTX_PFX, id_tag)
            obj_id = context.get(id_tag)
        else:
            obj_id = context.get(id_tag)
            context.update({tag4kids: obj_id})  # store for kids

        value = context.get("_LEAFVALUE")
        if type(value) == list and len(value) == 1:
            value = value[0]
        if type(value) != list and _a["key_is_list"] == True:
            # wrap back
            value = [value]
        if type(value) == list and _a["key_is_list"] == False:
            # force scalar
            value = ",".join(value)
        context_data = context.data.copy()
        context_data["_LEAFVALUE"] = value

        itag = None
        data = {}
        if _a["key"]:
            data = {_a["key"]: value}
        elif value is not None:
            itag = "_TECHVAL4IGNORE_"
            data[itag] = value
        if _a["addon"]:
            data.update(_a["addon"])

        data, _ = self.interpolate(data, context_data)
        if data is None:
            return

        if itag is not None:
            if itag in data:
                del data[itag]
            else:
                return

        context.global_context.add(obj_tag, obj_id, _a["path"], data, force=self._OUTPUT_FORCE_SUB)


class JsonSubRule(JsonRule):
    NAME = "JSON_SUB"
    _RULES = JsonRule.RulesType()
    _OUTPUT_FORCE_SUB = True
