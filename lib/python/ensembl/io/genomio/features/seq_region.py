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

from dataclasses import asdict, dataclass, field
from typing import Dict, List


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
    name: str
    length: int
    seq_region_id: int = 0
    coord_system_level: str = ""
    synonyms: List[SeqRegionSynonym] = field(default_factory=list)
    attributes: List[SeqRegionAttribute] = field(default_factory=list)
    karyotype_bands: List[KaryotypeBand] = field(default_factory=list)

    def get_seq_region_id(self):
        return self.seq_region_id
    
    def to_brc_dict(self) -> Dict:
        seqr_dict = {
            "name": self.name,
            "length": self.length,
            "coord_system_level": self.coord_system_level,
            "synonyms": [asdict(syn) for syn in self.synonyms],
            #"karyotype": asdict(band) for band in self.karyotype_bands)
        }
        attrib_dict = {attrib.code: attrib.value for attrib in self.attributes}
        
        attribs_to_add = {
            'BRC4_seq_region_name': 'BRC4_seq_region_name',
            'EBI_seq_region_name': 'EBI_seq_region_name',
            'coord_system_tag': 'coord_system_level',
            'sequence_location': 'location',
            'codon_table': 'codon_table',
            'circular_serq': 'circular',
        }
        for key, name in attribs_to_add.items():
            if key in attrib_dict:
                value = attrib_dict[key]
                seqr_dict[name] = value
        
        return seqr_dict