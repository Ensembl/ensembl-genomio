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
"""Construct a seq_region metadata file from INSDC files."""

__all__ = [
    "SeqRecordData",
    "SeqRegion",
    "UnknownMetadata",
]

import csv
from dataclasses import dataclass
import logging
from os import PathLike
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests

from ensembl.io.genomio.utils import get_json, print_json
from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

##############################################
SeqRegion = Dict[str, Any]

class SeqCollection(dict):
    mock: bool = False

    def from_gbff(self, gbff_path: Path | None = None) -> None:
        """Store seq_regions extracted from a GBFF file.
        If a seq_region with the same id exists in the collection, it will be replaced.
        """
        if gbff_file is None:
            return
        with open_gz_file(gbff_path) as gbff_file:
            for record in SeqIO.parse(gbff_file, "genbank"):
                record_data = SeqRecordData(record)
                seqr: SeqRegion = {"length": len(record.seq)}

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

                # Store the seq_region
                self[record.id] = seqr

    def from_report(self, report_path: Path | None = None, is_refseq: bool = False) -> None:
        """Store seq_regions extracted from an INSDC assembly report file.
        If a seq_region with the same id exists in the collection, it will be replaced.

        Args:
            report_path: Path to the sequence regions report file.
            is_refseq: True if the source of the report is RefSeq, false if INSDC.

        """
        if report_path is None:
            return
        # Get the report in a CSV format
        report_csv = self._report_to_csv(report_path)[0]
        # Feed the csv string to the CSV reader
        reader = csv.DictReader(report_csv.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE)
        # Create the seq_regions
        for row in reader:
            seq_region = self._make_seq_region(row, is_refseq)
            if seq_region:
                name = seq_region["name"]
                self[name] = seq_region

    @staticmethod
    def _report_to_csv(report_path: PathLike) -> Tuple[str, Dict]:
        """Returns an assembly report as a CSV string.

        Args:
            report_path: path to a seq_region file from INSDC/RefSeq

        Returns:
            The data as a string in CSV format, and the head metadata as a dictionary.

        """
        with open_gz_file(report_path) as report:
            data = ""
            metadata = {}
            last_head = ""
            for line in report:
                # Ignore header
                if line.startswith("#"):
                    # Get metadata values if possible
                    match = re.search("# (.+?): (.+?)$", line)
                    if match:
                        metadata[match.group(1)] = match.group(2)
                    last_head = line
                else:
                    if last_head:
                        data += last_head[2:].strip() + "\n"
                        last_head = ""
                    data += line
            return data, metadata

    @staticmethod
    def _make_seq_region(
        data: Dict,
        is_refseq: bool,
        synonym_map: Optional[Dict[str, str]] = None,
        molecule_location: Optional[Dict[str, str]] = None,
    ) -> SeqRegion:
        """Returns a sequence region from the information provided.

        An empty sequence region will be returned if not accession information is found.

        Args:
            data: a dict from the report representing one line, where the key is the column name.
            is_refseq: True if the source is RefSeq, false if INSDC.
            synonym_map: Map of INSDC report column names to sequence region field names.
            molecule_location: Map of sequence type to SO location.

        Raises:
            KeyError: If the sequence location is not recognised.
            UnknownMetadata: if the sequence role is not recognised.

        """
        if synonym_map is None:
            synonym_map = _SYNONYM_MAP
        if molecule_location is None:
            molecule_location = _MOLECULE_LOCATION
        seq_region = {}
        # Set accession as the sequence region name
        src = "RefSeq" if is_refseq else "GenBank"
        accession_id = data.get(f"{src}-Accn", "")
        if accession_id and (accession_id != "na"):
            seq_region["name"] = accession_id
        else:
            logging.warning(f'No {src} accession ID found for {data["Sequence-Name"]}')
            return {}
        # Add synonyms
        synonyms = []
        for field, source in synonym_map.items():
            if (field in data) and (data[field].casefold() != "na"):
                synonym = {"source": source, "name": data[field]}
                synonyms.append(synonym)
        if synonyms:
            synonyms.sort(key=lambda x: x["source"])
            seq_region["synonyms"] = synonyms
        # Add sequence length
        field = "Sequence-Length"
        if (field in data) and (data[field].casefold() != "na"):
            seq_region["length"] = int(data[field])
        # Add coordinate system and location
        seq_role = data["Sequence-Role"]
        # Scaffold?
        if seq_role in ("unplaced-scaffold", "unlocalized-scaffold"):
            seq_region["coord_system_level"] = "scaffold"
        # Chromosome? Check location
        elif seq_role == "assembled-molecule":
            seq_region["coord_system_level"] = "chromosome"
            location = data["Assigned-Molecule-Location/Type"].lower()
            # Get location metadata
            try:
                seq_region["location"] = molecule_location[location]
            except KeyError as exc:
                raise UnknownMetadata(f"Unrecognized sequence location: {location}") from exc
        else:
            raise UnknownMetadata(f"Unrecognized sequence role: {seq_role}")
        return seq_region

    def remove(self, to_exclude: list[str]):
        """Remove seq_regions based on a provided list of accessions."""
        for seq_name in to_exclude:
            if seq_name in self:
                del self[seq_name]
            else:
                logging.info(f'Cannot exclude seq not found: {seq_name}')
    
    def add_translation_table(
        self, location_codon: Optional[Dict[str, int]] = None
    ) -> None:
        """Adds the translation codon table to each sequence region (when missing) based on its location.

        Args:
            location_codon: Map of known codon tables for known locations.

        """
        if location_codon is None:
            location_codon = _LOCATION_CODON
        for seqr in self.values():
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
            return
        url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{str(taxon_id)}"
        response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
        response.raise_for_status()
        # In case we've been redirected, check for html opening tag
        if response.text.startswith("<"):
            raise ValueError(f"Response from {url} is not JSON")
        decoded = response.json()
        genetic_code = int(decoded.get("mitochondrialGeneticCode", 0))
        if genetic_code == 0:
            logging.warning(f"No mitochondria genetic code found for taxon {taxon_id}")
            return

        for seqr in self.values():
            if "codon_table" in seqr:
                continue
            if seqr.get("location", "") == "mitochondrial_chromosome":
                seqr["codon_table"] = genetic_code




##############################################


_SYNONYM_MAP = {
    "Assigned-Molecule": "INSDC",
    "GenBank-Accn": "GenBank",
    "RefSeq-Accn": "RefSeq",
    "Sequence-Name": "INSDC_submitted_name",
}
_MOLECULE_LOCATION = {
    "apicoplast": "apicoplast_chromosome",
    "chromosome": "nuclear_chromosome",
    "kinetoplast": "kinetoplast_chromosome",
    "linkage group": "linkage_group",
    "mitochondrion": "mitochondrial_chromosome",
    "plasmid": "plasmid",
}
_LOCATION_CODON = {"apicoplast_chromosome": 4}


class UnknownMetadata(Exception):
    """If a metadata if not supported or recognized."""


@dataclass
class SeqRecordData:
    record: SeqRecord

    def get_genbank_id(self) -> str | None:
        """Returns the GenBank accession from a given sequence record (if present).

        Only useful for RefSeq sequence records, where the GenBank accession is stored in a comment.

        Args:
            record: Sequence record.

        """
        genbank_id = None
        if "comment" in self.annotations:
            comment = str(self.annotations["comment"])
            comment = re.sub(r"[ \n\r]+", " ", comment)
            match = re.search(r"The reference sequence was derived from ([^\.]+)\.", comment)
            if match:
                genbank_id = match.group(1)
        return genbank_id


    def get_codon_table(self) -> Optional[int]:
        """Returns the codon table number from a given a GenBank sequence record (if present)."""
        table_number = None
        for feat in self.features:
            if "transl_table" in feat.qualifiers:
                table_number = int(feat.qualifiers["transl_table"][0])
                break
        return table_number


    def get_organelle(self, molecule_location: Optional[Dict] = None) -> Optional[str]:
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
        for feat in self.features:
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
        return self.annotations.get("topology", "") == "circular"


def prepare_seq_region_metadata(
    genome_file: PathLike,
    report_file: PathLike,
    dst_file: PathLike,
    gbff_file: Optional[PathLike] = None,
    to_exclude: Optional[List[str]] = None,
    mock_run: bool = False,
) -> None:
    """Prepares the sequence region metadata found in the INSDC/RefSeq report and GBFF files.

    The sequence region information is loaded from both sources and combined. Elements are added/excluded
    as requested, and the final sequence region metadata is dumped in a JSON file that follows the schema
    defined in "src/python/ensembl/io/genomio/data/schemas/seq_region.json".

    Args:
        genome_file: Genome metadata JSON file path.
        report_file: INSDC/RefSeq sequences report file path to parse.
        gbff_file: INSDC/RefSeq GBFF file path to parse.
        dst_file: JSON file output for the processed sequence regions JSON.
        to_exclude: Sequence region names to exclude.
        mock_run: Do not call external taxonomy service.

    """
    genome_data = get_json(genome_file)
    dst_file = Path(dst_file)
    is_refseq = genome_data["assembly"]["accession"].startswith("GCF_")

    seqs = SeqCollection(mock=mock_run)
    seqs.from_report(report_file, is_refseq)
    seqs.from_gbff(gbff_file)

    # Exclude seq_regions from a list
    if to_exclude is not None:
        seqs.remove(to_exclude)

    # Add translation and mitochondrial codon tables
    seqs.add_translation_table()
    seqs.add_mitochondrial_codon_table(genome_data["species"]["taxonomy_id"])

    # Print out the file
    print_json(dst_file, seqs)


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description="Construct a sequence region metadata file from INSDC files.")
    parser.add_argument_src_path("--genome_file", required=True, help="Genome metadata JSON file")
    parser.add_argument_src_path(
        "--report_file", required=True, help="INSDC/RefSeq sequences report file to parse"
    )
    parser.add_argument_src_path("--gbff_file", help="INSDC/RefSeq GBFF file to parse")
    parser.add_argument_dst_path(
        "--dst_file", default="seq_region.json", help="Output JSON file for the processed sequence regions"
    )
    parser.add_argument(
        "--to_exclude", nargs="*", metavar="SEQ_REGION_NAME", help="Sequence region names to exclude"
    )
    parser.add_argument("--mock_run", action="store_true", help="Do not call external APIs")
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    prepare_seq_region_metadata(
        genome_file=args.genome_file,
        report_file=args.report_file,
        dst_file=args.dst_file,
        gbff_file=args.gbff_file,
        to_exclude=args.to_exclude,
        mock_run=args.mock_run,
    )
