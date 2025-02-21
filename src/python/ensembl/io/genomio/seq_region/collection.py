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
"""SeqCollection as a collection of seq_regions metadata."""

__all__ = [
    "SeqCollection",
    "SeqRegionDict",
]

import logging
from pathlib import Path
from typing import Any, Mapping, TypeAlias

from Bio import SeqIO
import requests

from ensembl.io.genomio.seq_region.gbff import GBFFRecord
from ensembl.io.genomio.seq_region.exceptions import UnknownMetadata
from ensembl.io.genomio.seq_region.mappings import SYNONYM_MAP, MOLECULE_LOCATION, LOCATION_CODON
from ensembl.io.genomio.seq_region.report import ReportRecord
from ensembl.utils.archive import open_gz_file

SeqRegionDict: TypeAlias = dict[str, Any]


class SeqCollection:
    """Represent a collection of seq_regions metadata."""

    mock: bool
    seqs: dict

    def __init__(self, mock: bool = False) -> None:
        self.seqs = {}
        self.mock = mock

    def from_gbff(self, gbff_path: Path) -> None:
        """Store seq_regions extracted from a GBFF file.

        If a seq_region with the same ID exists in the collection, it will be replaced.
        """
        with open_gz_file(gbff_path) as gbff_file:
            for record in SeqIO.parse(gbff_file, "genbank"):
                record_data = GBFFRecord(record)
                new_seq: SeqRegionDict = self.make_seqregion_from_gbff(record_data)
                name = record.id
                merged_seq = self._merge(new_seq, self.seqs.get(name, {}))
                self.seqs[name] = merged_seq

    def _merge(self, source: SeqRegionDict, destination: SeqRegionDict) -> SeqRegionDict:
        """Merge a source dict in a destination dict (only add values or append to lists)."""
        if not destination:
            return source
        for key, value in source.items():
            if isinstance(value, list):
                destination[key] += value
            else:
                destination[key] = value

        return destination

    @staticmethod
    def make_seqregion_from_gbff(record_data: GBFFRecord) -> SeqRegionDict:
        """Returns a seq_region dict extracted from a GBFF record."""
        seqr: SeqRegionDict = {"length": len(record_data.record.seq)}  # type: ignore[arg-type]

        if record_data.is_circular():
            seqr["circular"] = True

        # Is there a genetic code defined?
        codon_table = record_data.get_codon_table()
        if codon_table is not None:
            seqr["codon_table"] = codon_table

        # Is it an organelle?
        location = record_data.get_organelle()
        if location is not None:
            seqr["location"] = location

        # Is there a comment stating the Genbank record this is based on?
        genbank_id = record_data.get_genbank_id()
        if genbank_id is not None:
            seqr["synonyms"] = [{"source": "INSDC", "name": genbank_id}]

        return seqr

    def from_report(self, report_path: Path, is_refseq: bool = False) -> None:
        """Store seq_regions extracted from an INSDC assembly report file.

        If a seq_region with the same id exists in the collection, it will be replaced.

        Args:
            report_path: Path to the sequence regions report file.
            is_refseq: True if the source of the report is RefSeq, false if INSDC.

        """
        report = ReportRecord(report_path)
        for seq_data in report.reader:
            new_seq = self.make_seq_region_from_report(seq_data, is_refseq)
            name = new_seq["name"]
            merged_seq = self._merge(new_seq, self.seqs.get(name, {}))
            self.seqs[name] = merged_seq

    @staticmethod
    def make_seq_region_from_report(
        seq_data: dict[str, Any],
        is_refseq: bool,
        synonym_map: Mapping[str, str] = SYNONYM_MAP,
        molecule_location: Mapping[str, str] = MOLECULE_LOCATION,
    ) -> SeqRegionDict:
        """Returns a sequence region from the information provided.

        An empty sequence region will be returned if no accession information is found.

        Args:
            data: Dict from the report representing one line, where the key is the column name.
            is_refseq: True if the source is RefSeq, false if INSDC.
            synonym_map: Map of INSDC report column names to sequence region field names.
            molecule_location: Map of sequence type to SO location.

        Raises:
            UnknownMetadata: If the sequence role or location is not recognised.

        """
        seq_region = {}

        # Set accession as the sequence region name
        src = "RefSeq" if is_refseq else "GenBank"
        accession_id = seq_data.get(f"{src}-Accn", "")
        if not accession_id or (accession_id == "na"):
            logging.warning(f'No {src} accession ID found for {seq_data["Sequence-Name"]}')
            return {}
        seq_region["name"] = accession_id

        # Add synonyms
        synonyms = []
        for field, source in synonym_map.items():
            if (field in seq_data) and (seq_data[field].casefold() != "na"):
                synonym = {"source": source, "name": seq_data[field]}
                synonyms.append(synonym)
        synonyms.sort(key=lambda x: x["source"])
        seq_region["synonyms"] = synonyms

        # Add sequence length
        field = "Sequence-Length"
        if (field in seq_data) and (seq_data[field].casefold() != "na"):
            seq_region["length"] = int(seq_data[field])

        # Add coordinate system and location
        seq_role = seq_data["Sequence-Role"]
        # Scaffold?
        if seq_role in ("unplaced-scaffold", "unlocalized-scaffold"):
            seq_region["coord_system_level"] = "scaffold"
        # Chromosome? Check location
        elif seq_role == "assembled-molecule":
            seq_region["coord_system_level"] = "chromosome"
            location = seq_data["Assigned-Molecule-Location/Type"].lower()
            # Get location metadata
            try:
                seq_region["location"] = molecule_location[location]
            except KeyError as exc:
                raise UnknownMetadata(f"Unrecognized sequence location: {location}") from exc
        else:
            raise UnknownMetadata(f"Unrecognized sequence role: {seq_role}")

        return seq_region

    def remove(self, to_exclude: list[str]) -> None:
        """Remove seq_regions based on a provided list of accessions."""
        for seq_name in to_exclude:
            if seq_name in self.seqs:
                del self.seqs[seq_name]
            else:
                logging.info(f"Cannot exclude seq not found: {seq_name}")

    def add_translation_table(self, location_codon: Mapping[str, int] = LOCATION_CODON) -> None:
        """Adds the translation codon table to each sequence region (when missing) based on its location.

        Args:
            location_codon: Map of known codon tables for known locations.

        """
        for seqr in self.seqs.values():
            if "codon_table" in seqr:
                continue
            if seqr.get("location", "") in location_codon:
                seqr["codon_table"] = location_codon[seqr["location"]]

    def add_mitochondrial_codon_table(self, taxon_id: int) -> None:
        """Adds the mitochondrial codon table to each sequence region (when missing) based on the taxon ID.

        If no mitochondrial genetic code can be found for the given taxon ID nothing will be changed.

        Args:
            taxon_id: The species taxon ID.

        """
        if self.mock:
            logging.info("Skip mitochondrial codon table: mock")
            return
        if not taxon_id:
            logging.info("Skip mitochondrial codon table: no taxon_id to use")
            return

        url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{str(taxon_id)}"
        response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
        response.raise_for_status()
        # In case we have been redirected, check for HTML opening tag
        if response.text.startswith("<"):
            raise ValueError(f"Response from {url} is not JSON")
        decoded = response.json()
        genetic_code = int(decoded.get("mitochondrialGeneticCode", 0))
        if genetic_code == 0:
            logging.warning(f"No mitochondria genetic code found for taxon {taxon_id}")
            return

        for seqr in self.seqs.values():
            if ("codon_table" not in seqr) and (seqr.get("location", "") == "mitochondrial_chromosome"):
                seqr["codon_table"] = genetic_code

    def to_list(self) -> list[SeqRegionDict]:
        """Returns the sequences as a simple list of `SeqRegionDict` objects."""
        return list(self.seqs.values())
