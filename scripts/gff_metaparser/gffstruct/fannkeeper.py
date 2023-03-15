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


import json
import sys

from collections import defaultdict, OrderedDict

from .basekeeper import BaseKeeper


class FannKeeper(BaseKeeper):
    # storing result functional annotation object
    def __init__(self):
        self._data = defaultdict(lambda: defaultdict(OrderedDict))

    def get(self, obj_tag, obj_id, path):
        if not obj_tag:
            return
        if obj_id and type(obj_id) == list:
            obj_id = obj_id[0]
        if not obj_id:
            return
        # find a place to get_from
        top = self._data[obj_tag][obj_id]
        if not path:
            return top
        for p in path.split("/"):
            if p not in top:
                return None
            top = top[p]
        return top

    def add(self, obj_tag, obj_id, path, value, force=False):
        if not value:
            return

        if obj_id and type(obj_id) == list:
            obj_id = obj_id[0]

        if not obj_id:
            return

        # find a place to insert to
        top = self._data[obj_tag][obj_id]
        if path:
            for p in path.split("/"):
                if type(p) == list:
                    Except("To many list levels for %s:%s/%s" % (obj_tag, obj_id, path))
                p_is_list = False
                if p.startswith("@"):
                    p_is_list = True
                    p = p[1:]
                if p not in top:
                    if p_is_list:
                        top[p] = list()
                    else:
                        top[p] = dict()
                top = top[p]

        # update top
        if type(top) == list:
            if not force:
                if not self.list_has(top, value):
                    top.append(value)
            else:
                top.clear()
                top.append(value)
            return
        for k, v in value.items():
            if not v:
                continue
            if k not in top:
                top[k] = v
                continue
            if type(v) == list:
                if not force:
                    top[k] += v
                    top[k] = self.uniq_list(top[k])
                else:
                    top[k] = [v]
                continue
            if not force:
                continue
            top[k] = v
        return

    def uniq_list(self, lst):
        # slow
        if len(lst) <= 1:
            return lst
        return [json.loads(p) for p in frozenset(map(lambda x: json.dumps(x, sort_keys=True), lst))]

    def list_has(self, lst, obj):
        # slow
        if not lst or obj is None or len(lst) < 1:
            return False
        prj = json.dumps(obj, sort_keys=True)
        return prj in frozenset(map(lambda x: json.dumps(x, sort_keys=True), lst))

    def dump_json(self, out_file, maps=None, dump_filter=None):
        if not out_file:
            return
        vals = []
        for tag in sorted(self._data.keys()):
            no_stashed = []
            for h in self._data[tag].values():
                x = h
                if "_STASH" in x:
                    x = h.copy()
                    x.pop("_STASH")
                no_stashed.append(x)
            vals += list(filter(dump_filter, no_stashed))
        json.dump(vals, out_file, indent=2, sort_keys=True)

    def dump(self, out_file, maps=None, dump_filter=None):
        self.dump_json(out_file, maps, dump_filter)
