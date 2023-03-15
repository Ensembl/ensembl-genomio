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


# finding most commomn prefix
class PfxTr:
    def __init__(self):
        self.letters = dict()
        self.count = int(0)

    def add(self, s):
        if s is None or s == "":
            return
        root = self
        for c in [""] + list(s):
            if c not in root.letters:
                root.letters[c] = PfxTr()
            root = root.letters[c]
            root.count += 1

    def get_max(self):
        (pfx, max_cnt) = ("", 0)
        root = self
        while root.letters:
            c = max(root.letters.keys(), key=lambda x: root.letters[x].count)
            if root.letters[c].count < max_cnt:
                return (pfx, max_cnt)
            max_cnt = root.letters[c].count
            pfx += c
            root = root.letters[c]
        return (pfx, max_cnt)
