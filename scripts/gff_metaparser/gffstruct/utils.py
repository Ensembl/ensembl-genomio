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


import gzip
import re
import sys

from collections import defaultdict
from Bio import SeqIO


class IdTrimmer:
    def __init__(self, re_rules=dict()):
        self._rules = None
        self.compile_rules(re_rules)

    def compile_rules(self, re_rules):
        if re_rules:
            self._rules = {k.lower(): re.compile(v) for k, v in re_rules.items()}

    def supporting(self, type=None):
        if type is None or not self._rules:
            return False
        if "any" in self._rules:
            return True
        return type.lower() in self._rules

    def normalize(self, id_str, type=None, context=None):
        id_str = str(id_str)
        _type = type.lower()
        if not self._rules:
            return id_str
        if type is None:
            _type = "any"
        if _type in self._rules:
            return self._rules[_type].sub("", id_str)
        if "any" in self._rules:
            return self._rules["any"].sub("", id_str)
        return id_str

    def __call__(self, *args, **kwargs):
        return self.normalize(*args, **kwargs)


class PfxTrimmer(IdTrimmer):
    def __init__(self, trim_str):
        self._rules = None
        if not trim_str:
            return
        rules = defaultdict(list)
        for trim_pair in filter(None, trim_str.split(",")):
            tag, *pat = trim_pair.split(":", 1)
            if not pat:
                tag, pat = "ANY", [tag]
            pat = pat[0]
            if tag.endswith("!"):
                tag = tag[:-1]
            else:
                pat = re.escape(pat)
            rules[tag].append(pat)
        rules = {k: r"^(?:%s)" % ("|".join(rules.get("ANY", []) + v)) for k, v in rules.items()}
        # print(rules, file=sys.stderr)
        self.compile_rules(rules)


class ExtMapper:
    def __init__(self, tag, map_file=None, map_str=None):
        pass

    def map(self, val):
        pass

    def __call__(self, *args, **kwargs):
        return self.map(*args, **kwargs)


class SeqLenDict:
    def __init__(self, fna_file=None):
        self._len = None
        self.load_from_file(fna_file)

    def load_from_file(self, fasta_file):
        if not fasta_file:
            return
        _open = fasta_file.endswith(".gz") and gzip.open or open
        with _open(fasta_file, "rt") as fasta:
            fasta_parser = SeqIO.parse(fasta, "fasta")
            self._len = dict()
            for rec in fasta_parser:
                self._len[rec.name] = len(rec)

    def get_len(self, srid):
        if self._len and srid in self._len:
            return self._len[srid]
        return None

    def __call__(self, srid):
        return self.get_len(srid)


class UpdatingLen:
    def __init__(self, val, force_update=True):
        self._force_update = force_update
        self._update = False
        self._val = val
        if not self._val:
            self._update = True

    def is_updating(self):
        return self._update

    def update(self, val, stop_on_success=False):
        if not self._update:
            return
        if not val:
            return
        if not self._val or (self._val < val and self._force_update):
            self._val = val
        if self._val and stop_on_success:
            self._update = False

    def __call__(self, *args, **kwargs):
        return self.update(*args, **kwargs)

    def __int__(self):
        return self._val or int(0)

    def __len__(self):
        return self.__int__()

    def __str__(self):
        return str(self.__len__())
