#!env python3
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

from pathlib import Path
import re
import csv
from typing import Any, Dict, Tuple

from ensembl.utils.archive import open_gz_file


class SeqregionParser:
    """Parser of a seq_region report from INSDC/RefSeq.

    The main method of the Parser is get_report_regions, which returns a Dict of seq_regions,
    where the keys are the names.
    """

    synonym_map = {
        "Sequence-Name": "INSDC_submitted_name",
        "GenBank-Accn": "INSDC",
        "RefSeq-Accn": "RefSeq",
        "Assigned-Molecule": "GenBank",
    }

    molecule_location = {
        "chromosome": "nuclear_chromosome",
        "mitochondrion": "mitochondrial_chromosome",
        "apicoplast": "apicoplast_chromosome",
        "plasmid": "plasmid",
        "linkage group": "linkage_group",
        "kinetoplast": "kinetoplast",
    }

    def get_report_regions(
        self, report_path: Path, accession: str, use_refseq: bool = False
    ) -> Dict[str, dict]:
        """Get seq_region data from report file.

        Args:
            report_path: Path to the INSDC seq_region report.
            use_refseq: Expect a RefSeq seq_region report.
        Returns:
            A dict of seq_regions dicts, with their name as the key
        """
        if accession.startswith("GCF"):
            use_refseq = True
        # Get the report in a CSV format, easier to manipulate
        report_csv, metadata = self.report_to_csv(report_path)

        # Feed the csv string to the CSV reader
        reader = csv.DictReader(report_csv.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE)

        # Metadata
        assembly_level = "contig"
        if "Assembly level" in metadata:
            assembly_level = metadata["Assembly level"].lower()

        # Create the seq_regions
        seq_regions = {}
        for row in reader:
            seq_region = self.make_seq_region(row, assembly_level, use_refseq)
            name = seq_region["name"]
            seq_regions[name] = seq_region

        return seq_regions

    def report_to_csv(self, report_path: Path) -> Tuple[str, dict]:
        """Load an assembly report as a csv string.

        Args:
            report_path: Path to the INSDC seq_region report.
        Returns:
            The csv as a string, and the head metadata as a dict.
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
                    continue
                else:
                    if last_head:
                        data += last_head[2:].strip() + "\n"
                        last_head = ""
                    data += line
            return data, metadata

    def make_seq_region(self, row: Dict[str, str], assembly_level: str, use_refseq: bool) -> Dict[str, Any]:
        """From a row of the report, create one seq_region dict.

        Args:
            row: A seq_region row from the INSDC report.
            assembly_level: what level is the seq_region (chromosome, contig, etc.)
            use_refseq: Expect a RefSeq seq_region report.

        Returns:
            A seq_region dict.
        """
        seq_region: Dict[str, Any] = {}

        # Map the fields to their synonym name
        synonym_map = self.synonym_map
        molecule_location = self.molecule_location

        # Synonyms
        synonyms = []
        for field, source in synonym_map.items():
            if field in row and row[field].lower() != "na":
                synonym = {"source": source, "name": row[field]}
                synonyms.append(synonym)
        if len(synonyms) > 0:
            synonyms.sort(key=lambda x: x["source"])
            seq_region["synonyms"] = synonyms

        # Length
        field = "Sequence-Length"
        if field in row and row[field].lower() != "na":
            seq_region["length"] = int(row[field])

        if use_refseq:
            if "RefSeq-Accn" in row:
                seq_region["name"] = row["RefSeq-Accn"]
            else:
                raise Exception(f"No RefSeq name for {row['Sequence-Name']}")

        else:
            if "GenBank-Accn" in row:
                seq_region["name"] = row["GenBank-Accn"]
            else:
                raise Exception(f"No INSDC name for {row['Sequence-Name']}")

        # Coord system and location
        seq_role = row["Sequence-Role"]

        # Scaffold?
        if seq_role in ("unplaced-scaffold", "unlocalized-scaffold", "alt-scaffold"):
            seq_region["coord_system_level"] = "scaffold"

        # Chromosome? Check location
        elif seq_role == "assembled-molecule":
            seq_region["coord_system_level"] = "chromosome"
            location = row["Assigned-Molecule-Location/Type"].lower()

            # Get location metadata
            if location in molecule_location:
                seq_location = molecule_location[location]
                seq_region["location"] = seq_location
            else:
                raise Exception(f"Unrecognized sequence location: {location}")
        else:
            raise Exception(f"Unrecognized sequence role: {seq_role}")
        return seq_region
