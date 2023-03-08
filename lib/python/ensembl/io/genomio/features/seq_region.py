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
from pathlib import Path
from typing import ClassVar, Dict, List, Set


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

    def __post_init__(self):
        self._update_source()

    def _update_source(self) -> None:
        """Using an external map, update the source name"""
        if self._db_map and self.source in self._db_map:
            self.source = self._db_map[self.source]

    @classmethod
    def set_map(cls, map_file: Path) -> None:
        """Class method, set up the map for all SeqRegion objects"""
        cls._db_map = dict()
        with map_file.open("r") as map_fh:
            for line in map_fh:
                line = line.rstrip()
                if line.startswith("#") or line.startswith(" ") or line == "":
                    continue
                parts = line.split("\t")
                if not parts[0] or not parts[1]:
                    raise Exception(f"External db file is not formatted correctly for: {line}")
                else:
                    cls._db_map[parts[1]] = parts[0]

    def to_brc_dict(self):
        syn_dict = {
            "name": self.synonym,
            "source": self.source
        }
        return syn_dict


@dataclass
class SeqRegionAttribute:
    value: str = ""
    code: str = ""
    integer_fields: ClassVar[Set] = {"codon_table"}
    bool_fields: ClassVar[Set] = {"circular_seq", "non_ref"}

    @classmethod
    def is_integer_field(cls, key: str) -> bool:
        return (key in cls.integer_fields)

    @classmethod
    def is_bool_field(cls, key: str) -> bool:
        return (key in cls.bool_fields)


@dataclass
class KaryotypeBand:
    band: str
    seq_region_start: int
    seq_region_end: int
    stain: str

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
class Provider:
    name: str
    url: str


@dataclass
class AddedSequence:
    accession: str
    assembly_provider: Provider
    annotation_provider: Provider


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

    def _get_added_sequence(self) -> AddedSequence:
        attrib_dict = self._get_attrib_dict()
        accession = attrib_dict.get("added_seq_accession")
        if accession:
            assembly_provider = Provider(
                name=attrib_dict.get("added_seq_asm_pr_nam"),
                url=attrib_dict.get("added_seq_asm_pr_url")
            )
            annotation_provider = Provider(
                name=attrib_dict.get("added_seq_ann_pr_nam"),
                url=attrib_dict.get("added_seq_ann_pr_url")
            )
            added_seq = AddedSequence(
                accession=accession,
                assembly_provider=assembly_provider,
                annotation_provider=annotation_provider
            )
            return added_seq

    def to_brc_dict(self) -> Dict:
        seqr_dict = {
            "name": self.name,
            "length": self.length,
            "coord_system_level": self._get_coord_system_level(),
        }
        if self.synonyms:
            seqr_dict["synonyms"] = [syn.to_brc_dict() for syn in self.synonyms]

        if self.karyotype_bands:
            seqr_dict["karyotype"] = [band.to_brc_dict() for band in self.karyotype_bands]

        added_seq = self._get_added_sequence()
        if added_seq:
            seqr_dict["added_sequence"] = asdict(added_seq)

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
                if SeqRegionAttribute.is_integer_field(key):
                    if value.isdigit():
                        value = int(value)
                    else:
                        raise ValueError(f"Integer expected for field '{key}', got: '{value}'")

                if SeqRegionAttribute.is_bool_field(key):
                    if value:
                        value = True
                    else:
                        value = False

                seqr_dict[name] = value

        return seqr_dict

    def is_top_level(self):
        for attrib in self.attributes:
            if attrib.code == "toplevel" and attrib.value == "1":
                return True
        return False
