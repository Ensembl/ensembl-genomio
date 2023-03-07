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
class CoordSystem:
    coord_system_id: int
    species_id: int
    name: str
    version: str
    attrib: List[str] = field(default_factory=list)


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
    band: str
    seq_region_start: int
    seq_region_end: int
    stain: str
    structure: str = ""

    def _get_structure(self):
        struct_dict = {
            "TEL": "telomere",
            "ACEN": "centromere"
        }
        return struct_dict.get(self.stain, "")

    def to_brc_dict(self):
        kar_dict = {
            "name": self.band,
            "start": self.seq_region_start,
            "end": self.seq_region_end
        }
        if self.stain:
            kar_dict["stain"] = self.stain
            kar_dict["structure"] = self._get_structure()

        return kar_dict


@dataclass
class SeqRegion:
    name: str
    length: int
    seq_region_id: int = 0
    coord_system_id: int = 0
    coord_system: CoordSystem = None
    synonyms: List[SeqRegionSynonym] = field(default_factory=list)
    attributes: List[SeqRegionAttribute] = field(default_factory=list)
    karyotype_bands: List[KaryotypeBand] = field(default_factory=list)

    def get_seq_region_id(self):
        return self.seq_region_id

    def _get_attrib_dict(self):
        return {attrib.code: attrib.value for attrib in self.attributes}

    def _get_coord_system_level(self):
        coord_level = self.coord_system.name
        if coord_level == 'primary_assembly':
            attrib = self._get_attrib_dict()
            coord_tag = attrib.get("coord_system_tag")
            if coord_tag:
                coord_level = coord_tag
            else:
                coord_level = ""

        return coord_level

    def to_brc_dict(self) -> Dict:
        seqr_dict = {
            "name": self.name,
            "length": self.length,
            "coord_system_level": self._get_coord_system_level(),
        }
        if self.synonyms:
            seqr_dict["synonyms"] = [asdict(syn) for syn in self.synonyms]

        if self.karyotype_bands:
            seqr_dict["karyotype"] = [band.to_brc_dict() for band in self.karyotype_bands]

        attrib_dict = self._get_attrib_dict()

        attribs_to_add = {
            'BRC4_seq_region_name': 'BRC4_seq_region_name',
            'EBI_seq_region_name': 'EBI_seq_region_name',
            'coord_system_tag': 'coord_system_level',
            'sequence_location': 'location',
            'codon_table': 'codon_table',
            'circular_seq': 'circular',
            'non_ref': 'non_ref',
        }
        for key, name in attribs_to_add.items():
            if key in attrib_dict:
                value = attrib_dict[key]
                if value.isdigit():
                    value = int(value)
                seqr_dict[name] = value

        return seqr_dict
