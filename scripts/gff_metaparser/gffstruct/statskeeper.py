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


import re
import sys

from collections import defaultdict, OrderedDict

from .basekeeper import BaseKeeper
from .prefixtree import PfxTr


class StatsKeeper(BaseKeeper):
    # storing stats
    def __init__(self, detailed=None, no_longest_pfx=False, no_id_stats_rules=[]):
        self._data = defaultdict(lambda: defaultdict(dict))
        self._detailed = detailed
        self._no_longest_pfx = no_longest_pfx
        self._no_id_stats_rules = frozenset(no_id_stats_rules)
        self._idstat = defaultdict(lambda: defaultdict(int))
        self._tr_d_0 = str.maketrans("123456789", "0" * 9)
        self._tr_00_0 = re.compile("0+")

    def add(self, rule_name, context):
        _fulltag = context.get("_FULLTAG")
        tagstat = self._data[rule_name][_fulltag]

        if "counts" not in tagstat:
            tagstat["counts"] = 0
        tagstat["counts"] += 1

        type_id_coords = context.get_to_root(getter=lambda x: self.type_id_coords_from_ctx(x))

        if self._detailed:
            print("rule %s for %s: %s " % (rule_name, _fulltag, type_id_coords), file=self._detailed)

        if not self._no_longest_pfx:
            if "pfx" not in tagstat:
                tagstat["pfx"] = defaultdict(PfxTr)
            for _type, _id, _ in type_id_coords:
                if _id and _id != ".":
                    tagstat["pfx"][_type].add(_id)
                    if not self._no_id_stats_rules or rule_name not in self._no_id_stats_rules:
                        self._idstat[_type][self.stem_id(_id)] += 1
                        self._idstat[_type]["_ALL"] += 1
        return

    def stem_id(self, s):
        pre = str(s).translate(self._tr_d_0)
        return self._tr_00_0.sub("0", pre)

    def summary(self, rule_name, out_file=None):
        if rule_name not in self._data:
            return None
        out = []
        for _tag in self._data[rule_name]:
            tagstat = self._data[rule_name][_tag]
            counts = str(tagstat.get("counts"))
            chain = str([{tp: pfx.get_max()} for tp, pfx in tagstat.get("pfx", {}).items()])
            out_line = "\t".join(["## stats", rule_name, _tag, counts, chain])
            if out_file:
                print(out_line, file=out_file)
            else:
                out.append(out_line)
        if out_file:
            return True
        return "\n".join(out)

    def dump(self, out_file, id_stats=True):
        for _rule in self._data:
            self.summary(_rule, out_file)
        if not id_stats:
            return
        for _type, _id_cnt in sorted(self._idstat.items(), key=lambda k: k[0]):
            for _id, _cnt in sorted(_id_cnt.items(), key=lambda k: -k[1]):
                if _id == "_ALL":
                    continue
                print("## stats\tID\t%s\t%s\t%s\t%s" % (_type, _id_cnt["_ALL"], _id, _cnt), file=out_file)

    def type_id_coords_from_ctx(self, context):
        if not context:
            return ".", ".", "."
        _type, _id, _seqid, _start, _end, _strand = map(
            lambda x: context.get("_" + x) or ".", "TYPE ID SEQID START END STRAND".split()
        )
        _strand = _strand == "." and "." or int(_strand) < 0 and "-" or "+"
        _coords = "%s:%s:%s-%s" % (_seqid, _strand, _start, _end)
        return _type, _id, _coords
