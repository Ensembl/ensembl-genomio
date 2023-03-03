#!/usr/bin/env python
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
"""Simple representation of a seq_region from an EnsEMBL database
"""

from dataclasses import dataclass, field
from typing import List


@dataclass
class SeqRegionSynonym:
    synonym: str = ""
    source: str = ""


@dataclass
class SeqRegionAttribute:
    value: str = ""
    code: str = ""


@dataclass
class KaryotypeBand:
    start: int
    end: int
    name: str
    stain: str
    structure: str


@dataclass
class SeqRegion:
    _seq_region_id: int
    name: str
    length: int
    coord_system_level: str = ""
    synonyms: List[SeqRegionSynonym] = field(default_factory=list)
    attributes: List[SeqRegionAttribute] = field(default_factory=list)
    karyotype_bands: List[KaryotypeBand] = field(default_factory=list)

    def get_seq_region_id(self):
        return self._seq_region_id
