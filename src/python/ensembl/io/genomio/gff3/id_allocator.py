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

from dataclasses import dataclass, field
import logging
import re
from typing import Dict, List, Optional, Set

from .features import GFFSeqFeature


class InvalidStableID(ValueError):
    """Raised when there is a problem with an stable ID."""


@dataclass
class StableIDAllocator:
    """Set of tools to check and allocate stable IDs."""

    # Multiple parameters to automate various fixes
    skip_gene_id_validation: bool = False
    min_id_length: int = 7
    current_id_number: int = 0
    make_missing_stable_ids: bool = True
    prefix: str = "TMP_"
    _loaded_ids: Set = field(default_factory=set)

    def set_prefix(self, genome: Dict) -> None:
        """Sets the ID prefix using the organism abbrev if it exists in the genome metadata."""
        try:
            org = genome["BRC4"]["organism_abbrev"]
        except KeyError:
            prefix = "TMP_PREFIX_"
        else:
            prefix = "TMP_" + org + "_"
        self.prefix = prefix

    def generate_gene_id(self) -> str:
        """Returns a new unique gene stable_id with a prefix.

        The ID is made up of a prefix and a number, which is auto incremented.

        """
        self.current_id_number += 1
        new_id = f"{self.prefix}{self.current_id_number}"
        return new_id

    def is_valid(self, stable_id: str) -> bool:
        """Checks that the format of a stable ID is valid.
        Args:
            stable_id: Stable ID to validate.
        """

        if self.skip_gene_id_validation:
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
        """Returns the stable ID after removing its prefix (if any).

        If more than one prefix may be found, only the first one is removed.

        Args:
            stable_id: Stable ID to process.
            prefixes: List of prefixes to search for.
        """

        for prefix in prefixes:
            if stable_id.startswith(prefix):
                return stable_id[len(prefix) :]
        return stable_id

    @staticmethod
    def generate_transcript_id(gene_id: str, number: int) -> str:
        """Returns a formatted transcript ID generated from a gene ID and number.
        Args:
            gene_id: Gene stable ID.
            number: Positive number.
        Raises:
            ValueError: If the number provided is not greater than zero.

        """
        if number < 1:
            raise ValueError("Number has to be a positive integer.")

        transcript_id = f"{gene_id}_t{number}"
        return transcript_id

    def normalize_cds_id(self, cds_id: str) -> str:
        """Returns a normalized version of the provided CDS ID.

        The normalisation implies to remove any unnecessary prefixes around the CDS ID. However, if
        the CDS ID is still not proper, an empty string will be returned.

        Args:
            cds_id: CDS ID to normalize.

        """

        prefixes = ["cds-", "cds:"]
        normalized_cds_id = StableIDAllocator.remove_prefix(cds_id, prefixes)

        # Special case: if the ID doesn't look like one, remove it - it needs to be regenerated
        if not self.is_valid(normalized_cds_id):
            return ""
        return normalized_cds_id

    def normalize_pseudogene_cds_id(self, pseudogene: GFFSeqFeature) -> None:
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

    def normalize_gene_id(self, gene: GFFSeqFeature, refseq: Optional[bool] = False) -> str:
        """Returns a normalized gene stable ID.

        Removes any unnecessary prefixes, but will generate a new stable ID if the normalized one is
        not recognized as valid.

        Args:
            gene: Gene feature to normalize.
        """
        prefixes = ["gene-", "gene:"]
        new_gene_id = StableIDAllocator.remove_prefix(gene.id, prefixes)

        is_valid = False
        # Special case for RefSeq: only valid Gene IDs are LOC*
        if refseq:
            if new_gene_id.startswith("LOC"):
                is_valid = True
        else:
            is_valid = self.is_valid(new_gene_id)

        if is_valid:
            return new_gene_id

        # In case the normalized gene ID is not valid, use the GeneID
        logging.debug(f"Gene ID is not valid: {new_gene_id}")
        qual = gene.qualifiers
        if "Dbxref" in qual:
            for xref in qual["Dbxref"]:
                (db, value) = xref.split(":")
                if db != "GeneID":
                    continue
                new_gene_id_base = f"{db}_{value}"
                new_gene_id = new_gene_id_base
                number = 1
                while new_gene_id in self._loaded_ids:
                    number += 1
                    new_gene_id = f"{new_gene_id_base}_{number}"
                    if number > 10:
                        raise InvalidStableID(f"Duplicate ID {new_gene_id_base} (up to {new_gene_id})")
                self._loaded_ids.add(new_gene_id)
                logging.debug(f"Using GeneID {new_gene_id} for stable_id instead of {gene.id}")
                return new_gene_id

        # Make a new stable_id
        if self.make_missing_stable_ids:
            new_gene_id = self.generate_gene_id()
            logging.debug(f"New ID: {new_gene_id} -> {new_gene_id}")
            return new_gene_id
        raise InvalidStableID(f"Can't use invalid gene id for {gene}")
