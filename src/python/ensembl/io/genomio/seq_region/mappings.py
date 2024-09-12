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
"""Seq region mappings."""

__all__ = [
    "SYNONYM_MAP",
    "MOLECULE_LOCATION",
    "LOCATION_CODON",
]

from types import MappingProxyType
from typing import Mapping


SYNONYM_MAP: Mapping[str, str] = MappingProxyType(
    {
        "Assigned-Molecule": "INSDC",
        "GenBank-Accn": "GenBank",
        "RefSeq-Accn": "RefSeq",
        "Sequence-Name": "INSDC_submitted_name",
    }
)
MOLECULE_LOCATION: Mapping[str, str] = MappingProxyType(
    {
        "apicoplast": "apicoplast_chromosome",
        "chromosome": "nuclear_chromosome",
        "kinetoplast": "kinetoplast_chromosome",
        "linkage group": "linkage_group",
        "mitochondrion": "mitochondrial_chromosome",
        "plasmid": "plasmid",
    }
)
LOCATION_CODON: Mapping[str, int] = MappingProxyType({"apicoplast_chromosome": 4})
