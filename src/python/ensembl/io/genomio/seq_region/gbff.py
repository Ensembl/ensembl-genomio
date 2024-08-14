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
"""A SeqRecord wrapper."""

__all__ = [
    "GBFFRecord",
]

from dataclasses import dataclass
import re

from Bio.SeqRecord import SeqRecord

from ensembl.io.genomio.seq_region.exceptions import UnknownMetadata

##############################################
_MOLECULE_LOCATION = {
    "apicoplast": "apicoplast_chromosome",
    "chromosome": "nuclear_chromosome",
    "kinetoplast": "kinetoplast_chromosome",
    "linkage group": "linkage_group",
    "mitochondrion": "mitochondrial_chromosome",
    "plasmid": "plasmid",
}


@dataclass
class GBFFRecord:
    """Wrapper around a SeqRecord object to extract specific data."""

    record: SeqRecord

    def get_genbank_id(self) -> str | None:
        """Returns the GenBank accession from a given sequence record (if present).

        Only useful for RefSeq sequence records, where the GenBank accession is stored in a comment.

        Args:
            record: Sequence record.

        """
        comment = str(self.record.annotations.get("comment", ""))
        if not comment:
            return
        comment = re.sub(r"[ \n\r]+", " ", comment)
        match = re.search(r"The reference sequence was derived from ([^\.]+)\.", comment)
        if not match:
            return
        return match.group(1)

    def get_codon_table(self) -> int | None:
        """Returns the codon table number from a given a GenBank sequence record (if present)."""
        table_number = None
        for feat in self.record.features:
            if "transl_table" in feat.qualifiers:
                table_number = int(feat.qualifiers["transl_table"][0])
                break
        return table_number

    def get_organelle(self, molecule_location: dict | None = None) -> str | None:
        """Returns the organelle location from the given GenBank record (if present).

        Args:
            record: GenBank sequence record.
            molecule_location: Map of sequence type to SO location.

        Raises:
            KeyError: If the location is not part of the controlled vocabulary.

        """
        if molecule_location is None:
            molecule_location = _MOLECULE_LOCATION
        location = None
        for feat in self.record.features:
            if "organelle" in feat.qualifiers:
                organelle = str(feat.qualifiers["organelle"][0])
                if not organelle:
                    break
                # Remove plastid prefix
                with_prefix = re.match(r"^(plastid|mitochondrion):(.+)$", organelle)
                if with_prefix:
                    organelle = with_prefix[2]
                # Get controlled name
                try:
                    location = molecule_location[organelle]
                except KeyError as exc:
                    raise UnknownMetadata(f"Unrecognized sequence location: {organelle}") from exc
                break
        return location

    def is_circular(self) -> bool:
        """Check if the record says that the sequence is circular."""
        return self.record.annotations.get("topology", "") == "circular"
