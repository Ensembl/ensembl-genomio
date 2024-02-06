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
"""Check and allocate IDs for gene features in a GFF3 file."""

__all__ = ["StableIDAllocator", "InvalidStableID"]

from dataclasses import dataclass
import logging
import re
from typing import List

from Bio.SeqFeature import SeqFeature


class InvalidStableID(ValueError):
    """Raised when there is a problem with an ID."""


@dataclass
class StableIDAllocator:
    """Set of tools to check and allocate stable IDs."""

    # Multiple parameters to automate various fixes
    validate_gene_id: bool = True
    min_id_length: int = 8
    current_id_number: int = 0
    make_missing_stable_ids: bool = True
    prefix = "TMP_"

    def generate_gene_id(self) -> str:
        """Returns a new unique gene stable_id with a prefix.

        The ID is made up of a prefix and a number, which is auto incremented.

        """
        self.current_id_number += 1
        new_id = f"{self.prefix}{self.current_id_number}"
        return new_id

    def is_valid(self, stable_id: str) -> bool:
        """Check that the format of a stable id is valid."""

        if not self.validate_gene_id:
            logging.debug(f"Validation deactivated by user: '{stable_id}' not checked")
            return True

        # Trna (from tRNAscan)
        if re.search(r"^Trna", stable_id, re.IGNORECASE):
            logging.debug(f"Stable ID is a tRNA from tRNA-scan: {stable_id}")
            return False

        # Coordinates
        if re.search(r"^.+:\d+..\d+", stable_id):
            logging.debug(f"Stable id is a coordinate: {stable_id}")
            return False

        # Special characters
        if re.search(r"[ |]", stable_id):
            logging.debug(f"Stable id contains special characters: {stable_id}")
            return False

        # Min length
        if len(stable_id) < self.min_id_length:
            logging.debug(f"Stable id is too short (<{self.min_id_length}) {stable_id}")
            return False

        return True

    @staticmethod
    def remove_prefix(stable_id: str, prefixes: List[str]) -> str:
        """Returns the identifier after removing all the prefixes found in it (if any)."""
        for prefix in prefixes:
            if stable_id.startswith(prefix):
                stable_id = stable_id[len(prefix) :]
                break
        return stable_id

    @staticmethod
    def generate_transcript_id(gene_id: str, number: int) -> str:
        """Use a gene ID and a number to make a formatted transcript ID."""

        transcript_id = f"{gene_id}_t{number}"
        return transcript_id

    def normalize_cds_id(self, cds_id: str) -> str:
        """Returns a normalised version of the provided CDS ID.

        The normalisation implies to remove any unnecessary prefixes around the CDS ID. However, if
        the CDS ID is still not proper, an empty string will be returned.

        """

        prefixes = ["cds-", "cds:"]
        cds_id = StableIDAllocator.remove_prefix(cds_id, prefixes)

        # Special case: if the ID doesn't look like one, remove it
        # It needs to be regenerated
        if not self.is_valid(cds_id):
            cds_id = ""

        return cds_id

    def normalize_pseudogene_cds_id(self, pseudogene: SeqFeature) -> None:
        """Normalizes every CDS ID of the provided pseudogene.

        Ensure each CDS from a pseudogene has a proper ID:
        - Different from the gene
        - Derived from the gene if it is not proper

        Args:
            pseudogene: Pseudogene feature.
        """
        for transcript in pseudogene.sub_features:
            for feat in transcript.sub_features:
                if feat.type == "CDS":
                    feat.id = self.normalize_cds_id(feat.id)
                    if feat.id in ("", pseudogene.id):
                        feat.id = f"{transcript.id}_cds"
                        feat.qualifiers["ID"] = feat.id

    def normalize_gene_id(self, gene: SeqFeature) -> str:
        """Returns a normalized gene stable ID.

        Removes any unnecessary prefixes, but will generate a new stable ID if the normalized one is
        not recognized as valid.

        Args:
            gene: Gene feature to normalize.
        """
        prefixes = ["gene-", "gene:"]
        new_gene_id = self.remove_prefix(gene.id, prefixes)

        # In case the normalized gene ID is not valid, use the GeneID
        if not self.is_valid(new_gene_id):
            logging.warning(f"Gene ID is not valid: {new_gene_id}")
            qual = gene.qualifiers
            if "Dbxref" in qual:
                for xref in qual["Dbxref"]:
                    (db, value) = xref.split(":")
                    if db == "GeneID":
                        new_gene_id = f"{db}_{value}"
                        logging.debug(f"Using GeneID {new_gene_id} for stable_id instead of {gene.id}")
                        return new_gene_id

            # Make a new stable_id
            if self.make_missing_stable_ids:
                new_id = self.generate_gene_id()
                logging.debug(f"New id: {new_gene_id} -> {new_id}")
                return new_id
            raise InvalidStableID(f"Can't use invalid gene id for {gene}")

        return new_gene_id
