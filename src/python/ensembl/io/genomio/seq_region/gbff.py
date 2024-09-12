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
"""A `SeqRecord` wrapper."""

__all__ = [
    "GBFFRecord",
]

from dataclasses import dataclass
import re
from typing import Mapping

from Bio.SeqRecord import SeqRecord

from ensembl.io.genomio.seq_region.mappings import MOLECULE_LOCATION
from ensembl.io.genomio.seq_region.exceptions import UnknownMetadata


@dataclass
class GBFFRecord:
    """Wrapper around a `SeqRecord` object to extract specific data."""

    record: SeqRecord

    def get_genbank_id(self) -> str | None:
        """Returns the GenBank accession from a given sequence record (if present).

        Only useful for RefSeq sequence records, where the GenBank accession is stored in a comment.

        Args:
            record: Sequence record.

        """
        comment = str(self.record.annotations.get("comment", ""))
        if not comment:
            return None
        comment = re.sub(r"[ \n\r]+", " ", comment)
        match = re.search(r"The reference sequence was derived from ([^\.]+)\.", comment)
        if not match:
            return None
        return match.group(1)

    def get_codon_table(self) -> int | None:
        """Returns the codon table number from a given a GenBank sequence record (if present)."""
        for feat in self.record.features:
            if "transl_table" in feat.qualifiers:
                return int(feat.qualifiers["transl_table"][0])
        return None

    def get_organelle(self, molecule_location: Mapping[str, str] = MOLECULE_LOCATION) -> str | None:
        """Returns the organelle location from the given GenBank record (if present).

        Args:
            record: GenBank sequence record.
            molecule_location: Map of sequence type to SO location.

        Raises:
            UnknownMetadata: If the location is not part of the controlled vocabulary.

        """
        location = None
        for feat in self.record.features:
            if "organelle" not in feat.qualifiers:
                continue
            organelle = str(feat.qualifiers["organelle"][0])
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
        """Returns True if the record says that the sequence is circular, False otherwise."""
        return self.record.annotations.get("topology", "") == "circular"
