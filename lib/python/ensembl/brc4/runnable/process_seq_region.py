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


import os
import re
from typing import Any, Dict, List, Tuple
import eHive
import gzip
import csv
import json

from Bio import SeqIO, SeqRecord
import requests

SeqRegion = Dict[str, Any]


class MissingDataError(Exception):
    """Used if some data is missing from the report file."""
    def __init__(self, report_path: str, accession: str, msg):
        report_msg = f"Can't get data for {accession} in report {report_path}"
        if msg:
            report_msg = f"{report_msg}: {msg}"
        self.msg = report_msg


class process_seq_region(eHive.BaseRunnable):
    """Runnable to load seq_regions metadata from INSDC/RefSeq reports and dump a json file

    The runnable does some checks, adds some things (like codon tables), and dumps the list of
    seq_regions in a json file that follows the schema defined in schema/seq_region_schema.json.

    Params:
        genome_data: Genome data following the schema/genome_schema.json.
        work_dir: Directory where the new file will be created.
        report: Path to the INSDC/RefSeq sequences report to parse.
        gbff: Path to the INSDC/RefSeq gbff file to parse.
        brc4_mode: Activate BRC4 mode (default).
        exclude_seq_regions: Array of seq_regions to not include in the new file.
    
    Predefined params:
        synonym_map: Map from the INSDC report column names to the seq_region field names.
        molecule_location: map the sequence type to an SO location.
        location_codon: map a location to a codon table.
    
    Dataflows:
        2: One dict with 2 keys:
            metadata_type: 'seq_region'.
            metadata_json: Path to the new seq_region json.
    
    Notes:
        BRC4 mode does the following:
            Add BRC4_seq_region_name and EBI_seq_region_name for each seq_region.
    """

    def param_defaults(self):
        return {
            "brc4_mode": True,
            "brc4_synonym_source": "GenBank",
            "brc4_synonym_source_alt": "RefSeq",
            "synonym_map": {
                "Sequence-Name": "INSDC_submitted_name",
                "GenBank-Accn": "GenBank",
                "RefSeq-Accn": "RefSeq",
                "Assigned-Molecule": "INSDC",
            },

            "molecule_location": {
                "chromosome": "nuclear_chromosome",
                "mitochondrion": "mitochondrial_chromosome",
                "apicoplast": "apicoplast_chromosome",
                "plasmid": "plasmid",
                "kinetoplast": "kinetoplast_chromosome",
                "linkage group": "linkage_group"
            },
            "location_codon": {
                "apicoplast_chromosome": 4
            },
            "exclude_seq_regions": [],
        }

    def run(self):
        genome_data = self.param('genome_data')
        work_dir = self.param('work_dir')
        report_path = self.param('report')
        gbff_path = self.param('gbff')
        brc4_mode = self.param('brc4_mode')

        # Create dedicated work dir
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final file name
        metadata_type = "seq_region"
        new_file_name = metadata_type + ".json"
        final_path = os.path.join(work_dir, new_file_name)
        
        use_refseq = self.param('accession').startswith("GCF_")

        # Get seq_regions data from report and gff3, and merge them
        report_regions = self.get_report_regions(report_path, use_refseq)
        gbff_regions = self.get_gbff_regions(gbff_path)
        seq_regions = self.merge_regions(report_regions, gbff_regions)

        # Exclude seq_regions from a list
        to_exclude = self.param("exclude_seq_regions")
        if to_exclude:
            seq_regions = self.exclude_seq_regions(seq_regions, to_exclude)
        
        # Setup the BRC4_seq_region_name
        if brc4_mode:
            seq_regions = self.add_brc4_ebi_name(seq_regions)

        # Guess translation table
        self.guess_translation_table(seq_regions)
        self.add_mitochondrial_codon_table(seq_regions, genome_data["species"]["taxonomy_id"])

        # Print out the file
        self.print_json(final_path, seq_regions)

        # Flow out the file and type
        output = {
            "metadata_type": metadata_type,
            "metadata_json": final_path
        }
        self.dataflow(output, 2)

    def exclude_seq_regions(self, seq_regions: List[SeqRegion],
                            to_exclude: List[str]) -> List[SeqRegion]:
        """Remove some seq_regions given as a list.
        
        Args:
            seq_regions: List of current seq_regions.
            to_exclude: List of seq_region names to exclude.
        
        Returns:
            A list of seq_regions with the ones from the exclusion list removed.
        """
        new_seq_regions = []
        for seqr in seq_regions:
            if "name" in seqr and seqr["name"] in to_exclude:
                print("Remove seq_region %s" % seqr["name"])
            else:
                new_seq_regions.append(seqr)
        return new_seq_regions

    def guess_translation_table(self, seq_regions: List[SeqRegion]) -> None:
        """Guess the codon table based on the location.

        Args:
            seq_regions: List of current seq_regions.
        
        Note:
            This relies on a map (location_codon) that gives the codon table,
            from a known location.
            The codon_table data is directly added to the SeqRegion if found.
        """
        location_codon = self.param("location_codon")
        
        for seqr in seq_regions:
            # Don't overwrite the existing codon table
            if "codon_table" in seqr:
                continue
            if "location" in seqr and seqr["location"] in location_codon:
                seqr["codon_table"] = location_codon[seqr["location"]]

    def add_brc4_ebi_name(self, seq_regions: List[SeqRegion]) -> List[SeqRegion]:
        """Use the INSDC seq_region name without version as the default BRC4 and EBI names.

        Args:
            seq_regions: List of current seq_regions.

        Note:
            The seq_region_names data are directly added to the SeqRegion if found.
        """
        source = self.param("brc4_synonym_source")
        source_alt = self.param("brc4_synonym_source_alt")
        
        new_seq_regions = []
        for seqr in seq_regions:
            main_name = ''
            alt_name = ''
            if "synonyms" in seqr:
                for syn in seqr["synonyms"]:
                    if syn["source"] == source:
                        main_name = syn["name"]
                    if syn["source"] == source_alt:
                        alt_name = syn["name"]
            
            name = ''
            if main_name:
                name = main_name
            elif alt_name:
                name = alt_name
            else:
                raise Exception(f"Can't set BRC4/EBI name for {seqr}")
            
            parts = name.partition(".")
            flat_name = parts[0]
            seqr["BRC4_seq_region_name"] = flat_name
            seqr["EBI_seq_region_name"] = flat_name
            new_seq_regions.append(seqr)

        return new_seq_regions 

    def add_mitochondrial_codon_table(self, seq_regions: List[SeqRegion], tax_id: int) -> None:
        """Add the codon table for mitochondria based on taxonomy.

        Args:
            seq_regions: List of current seq_regions.
            tax_id: The species taxonomy id.
        
        Note:
            The codon_table data is directly added to the SeqRegion if found.
        """
        tax_codon_table = self.get_mito_tax_codon_table(tax_id)
        if not tax_codon_table:
            return
        
        for seqr in seq_regions:
            if "location" in seqr and seqr["location"] == "mitochondrial_chromosome":
                if "codon_table" in seqr:
                    continue
                seqr["codon_table"] = tax_codon_table

    def get_mito_tax_codon_table(self, tax_id: int) -> int:
        """Get the codon table for mitochondria for a given taxon id.

        Args:
            tax_id: Taxonomy id.
        
        Returns:
            The codon table number if found. 0 otherwise.
        """
        server = "https://www.ebi.ac.uk"
        ext = "/ena/data/taxonomy/v1/taxon/tax-id/" + str(tax_id)
        genetic_code: int = 0

        r = requests.get(
            server + ext, headers={"Content-Type": "application/json"})
        decoded = r.json()

        try:
            if "mitochondrialGeneticCode" in decoded:
                genetic_code = int(decoded["mitochondrialGeneticCode"])
        except KeyError:
            print(f"No Mitochondria genetic code found for taxon {tax_id}")

        return genetic_code
    
    def print_json(self, path: str, data: Any) -> None:
        """Generic data json dumper to a file.
        
        Args:
            path: Path to the json to create.
            data: Any data to store.
         """
        with open(path, "w") as json_out:
            json_out.write(json.dumps(data, sort_keys=True, indent=4))
    
    def merge_regions(self,
                      regions1: Dict[str, SeqRegion],
                      regions2: Dict[str, SeqRegion]) -> List[SeqRegion]:
        """
        Merge seq_regions from different sources.

        Args:
            regions1: First Dict of SeqRegions, where the key is the name of the seq_region.
            regions2: Second Dict of SeqRegions.
        
        Returns:
            Alist of SeqRegions.
        """
        if not regions1:
            regions1 = {}
        if not regions2:
            regions2 = {}
        
        # Get all the names
        names1 = frozenset(regions1)
        names2 = frozenset(regions2)
        all_names = names1.union(names2)

        # Create the seq_regions, merge if needed
        seq_regions = []
        for name in all_names:
            
            seqr1 = {}
            seqr2 = {}
            if name in names1:
                seqr1 = regions1[name]
            if name in names2:
                seqr2 = regions2[name]

            final_seqr = {}
            if seqr1 and seqr2:
                final_seqr = self.merge_dicts(seqr1, seqr2)
            elif seqr1:
                final_seqr = seqr1
            elif seqr2:
                final_seqr = seqr2
            else:
                raise Exception(f"No seq_region found for {name}")

            seq_regions.append(final_seqr)

        seq_regions.sort(key=lambda x: x["name"])

        return seq_regions
    
    def merge_dicts(self, dict1: dict, dict2: dict) -> dict:
        """Merge 2 dicts. If they have the same key, the second dict overrides the value"""
        new_dict = {}

        for key in dict1:
            new_dict[key] = dict1[key]
        for key in dict2:
            new_dict[key] = dict2[key]
        
        return new_dict


    def get_gbff_regions(self, gbff_path: str) -> Dict[str, SeqRegion]:
        """Get seq_region data from the gbff file.

        Args:
            gbff_path: Gbff file path to use.
        
        Returns:
            A dict of SeqRegions, with their name as the key.
        """
        if not gbff_path:
            return {}
        
        seq_regions = {}
        _open = gbff_path.endswith(".gz") and gzip.open or open
        with _open(gbff_path, 'rt') as gbff_file:

            for record in SeqIO.parse(gbff_file, "genbank"):
                seqr: SeqRegion = {}

                seqr["length"] = len(record.seq)
                
                # Is the seq_region circular?
                annotations = record.annotations
                if "topology" in annotations and annotations["topology"] == "circular":
                    seqr["circular"] = True

                # Is there a genetic code defined?
                codon_table = self.get_codon_table(record)
                if codon_table:
                    seqr["codon_table"] = codon_table

                # Is it an organelle?
                location = self.get_organelle(record)
                if location:
                    seqr["location"] = location
                
                # Is there a comment stating the Genbank record this is based on?
                genbank_id = self.get_genbank_from_comment(record)
                if genbank_id:
                    seqr["synonyms"] = [
                        {"source": "INSDC", "name": genbank_id}
                    ]
                
                # Store the seq_region
                if seqr:
                    seq_regions[record.id] = seqr
        
        return seq_regions
    
    def get_genbank_from_comment(self, record: SeqRecord) -> str:
        """Given a genbank record, find a Genbank ID in the comment (if RefSeq).

        Args:
            record: The GenBank record to look into.
        
        Returns:
            A genbank accession as a string. Empty string if not found.
        """
        genbank = ""
        if "comment" in record.annotations:
            comment = record.annotations["comment"]
            comment = re.sub(r'[ \n\r]+', " ", comment)

            match = re.search(r'The reference sequence was derived from ([^\.]+)\.', comment)
            if match:
                genbank = match.group(1)
                
        return genbank
    
    def get_codon_table(self, record: SeqRecord) -> int:
        """Given a genbank record, seeks codon table features.

        Args:
            record: The GenBank record to look into.

        Returns:
            The codon table number if found, 0 otherwise.
        """

        table = 0

        for feat in record.features:
            if "transl_table" in feat.qualifiers:
                table = int(feat.qualifiers["transl_table"][0])
                break
        
        return table
    
    def get_organelle(self, record: SeqRecord) -> str:
        """Given a genbank record, look for an organelle field.

        Args:
            record: The GenBank record to look into.

        Returns:
            the organelle location (empty string if not found).
        """

        location = 0
        molecule_location = self.param("molecule_location")

        for feat in record.features:
            if "organelle" in feat.qualifiers:
                organelle = str(feat.qualifiers["organelle"][0])
                if not organelle:
                    break

                # Remove plastid prefix
                with_prefix = re.match(r'^(plastid|mitochondrion):(.+)$', organelle)
                if with_prefix:
                    organelle = with_prefix[2]

                # Get controlled name
                try:
                    location = molecule_location[organelle]
                except KeyError:
                    raise KeyError(f"Unrecognized sequence location: {organelle}")
                break
        
        return location

    def get_report_regions(self, report_path: str, use_refseq: bool) -> Dict[str, SeqRegion]:
        """Get seq_region data from report file.

        Args:
            report_path: Path to the seq_regions report from INSDC/RefSeq.
            use_refseq: True if using RefSeq, False if INSDC.
        
        Returns:
            A dict of SeqRegions, with their name as the key.
        """

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
            if not seq_region:
                continue
            name = seq_region["name"]
            seq_regions[name] = seq_region
        
        return seq_regions
    
    def make_seq_region(self, row: Dict, assembly_level: str, use_refseq: bool) -> SeqRegion:
        """From a row of the seq_region report, create one seq_region.

        Args:
            row: a dict from the report representing one line, where the key is the column name.
            assembly_level: what the whole assembly level is supposed to be (chromosome, scaffold)
            use_refseq: True if using RefSeq, False if INSDC.
        
        Returns:
            A single SeqRegion.
        """
        seq_region: SeqRegion = {}

        # Map the fields to their synonym name
        synonym_map = self.param("synonym_map")
        molecule_location = self.param("molecule_location")

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
        name = "length"
        if field in row and row[field].lower() != "na":
            seq_region[name] = int(row[field])
        
        refseq_id = ''
        gb_id = ''
        if "RefSeq-Accn" in row:
            refseq_id = row["RefSeq-Accn"]
            if refseq_id == "na":
                refseq_id = ''

        if "GenBank-Accn" in row:
            gb_id = row["GenBank-Accn"]
            if gb_id == "na":
                gb_id = ''

        if use_refseq:
            if refseq_id:
                seq_region["name"] = refseq_id
            else:
                print("No RefSeq name for %s" % row["Sequence-Name"])
                return {}
        elif gb_id:
            seq_region["name"] = gb_id
        else:
            print("No Genbank name for %s" % row["Sequence-Name"])
            return {}
        
        # Coord system and location
        seq_role = row["Sequence-Role"]
        
        # Scaffold?
        if seq_role in ("unplaced-scaffold", "unlocalized-scaffold"):
            seq_region["coord_system_level"] = "scaffold"
        
        # Chromosome? Check location
        elif seq_role == "assembled-molecule":
            seq_region["coord_system_level"] = "chromosome"
            location = row["Assigned-Molecule-Location/Type"].lower()
            
            # Get location metadata
            try:
                seq_region["location"] = molecule_location[location]
            except KeyError:
                raise KeyError(f"Unrecognized sequence location: {location}")
        else:
            raise Exception("Unrecognized sequence role: %s" % seq_role)
        return seq_region

    def report_to_csv(self, report_path: str) -> Tuple[str, dict]:
        """Load an assembly report as a csv string.

        Args:
            report_path: path to a seq_region file from INSDC/RefSeq
        
        Returns:
            The data as a string in csv format, and the head metadata as a dict.
        """

        _open = report_path.endswith(".gz") and gzip.open or open
        with _open(report_path, 'rt') as report:
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
