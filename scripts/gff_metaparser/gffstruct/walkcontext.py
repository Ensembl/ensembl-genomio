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


import copy
from collections import defaultdict


class WalkContext:
    def __init__(self, tag="", global_context=None, ctg_len_inferer=None):
        self.data = dict()
        self._tag = tag
        self.global_context = global_context
        self.ctg_len = ctg_len_inferer
        self.processed_rules = []
        self.prev = []
        self._top = None
        self._useful_leaves = []
        self.used_id_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))  # id depth count

    def snap(self):
        # shallow data copy
        self.prev.append(copy.copy(self.data))
        return self.prev[-1]

    def top(self, *feature):
        if len(feature) == 0:
            return self._top
        if feature[0]:
            self._top = feature[0]

    def tag(self, *val):
        if len(val) == 0:
            return self._tag
        self._tag = str(val[0])

    def update(self, *key_val, force_clean=False, **kwargs):
        # can have either dict as key or string with non-empty val
        if len(key_val) == 1 and isinstance(key_val[0], dict):
            for k, v in key_val[0].items():
                self.update(k, v, force_clean=force_clean)
        elif len(key_val) == 2:
            key, val = key_val
            if val is not None:
                self.data[key] = val
            elif force_clean and key in self.data:
                del self.data[key]
        # update from **kwargs
        for k, v in kwargs.items():
            self.update(k, v, force_clean=force_clean)

    def get(self, key, default=None):
        # global ???
        if key not in self.data:
            return default
        return self.data[key]

    def __getitem__(self, key):
        return self.get(key)

    def update_processed_rules(self, lst):
        self.processed_rules += lst

    def used_leaves(self, *leaves):
        if len(leaves) == 0:
            return copy.copy(self._useful_leaves)
        self._useful_leaves += leaves

    def add_fix(self, fix):
        pass

    def merge_fixes(self):
        pass

    def update_ctg_len(self, length):
        if self.ctg_len is not None:
            self.ctg_len(length)

    def get_to_root(self, getter=None, node=None):
        if not getter:
            return None
        res = []
        it = node
        if it is None:
            it = self.data
        while it:
            out = getter(it)
            if out:
                res.append(out)
            it = it.get("_PARENTCTX")
        return res[::-1]

    def run_to_root(self, updater=None, node=None):
        if not updater:
            return
        it = node
        if it is None:
            it = self.data
        while it:
            updater(it)
            it = it.get("_PARENTCTX")
        return

    def store_id_stats(self, _id, _type):
        self.used_id_stats[_id][0]["_ALL"] += 1
        self.used_id_stats[_id][0][_type.lower()] += 1
        # depth
        depth = self.data.get("_DEPTH", 0)
        self.used_id_stats[_id][depth]["_ALL"] += 1
        #
        self.used_id_stats[_id][depth][_type or "_ALL"] += 1

    def fix_id_sfx(self, raw_id, _type):
        if not raw_id:
            return raw_id
        if not raw_id in self.used_id_stats:
            self.store_id_stats(raw_id, _type)
            return raw_id
        # special cases
        depth = self.data.get("_DEPTH", 0)
        _type_lc = _type.lower()
        is_tr_type = "rna" in _type_lc or "transcript" in _type_lc
        # transcript
        if depth == 2 and is_tr_type:
            self.store_id_stats(raw_id, _type)
            self.store_id_stats(raw_id, "_TRANSCRIPT_LEVEL")
            _id = raw_id + "_t" + str(self.used_id_stats[raw_id][depth]["_TRANSCRIPT_LEVEL"])
            self.store_id_stats(_id, _type)
            return _id
        # cds -- multifeature
        if _type_lc == "cds":
            self.store_id_stats(raw_id, _type)
            allc = self.used_id_stats[raw_id][0]["_ALL"]
            typec = self.used_id_stats[raw_id][0][_type_lc]
            if typec == allc:
                self.store_id_stats(raw_id, _type)
                return raw_id
        # usual case
        self.store_id_stats(raw_id, _type)
        _id = raw_id + "_" + str(self.used_id_stats[raw_id][0]["_ALL"])
        self.store_id_stats(_id, _type)
        return _id
