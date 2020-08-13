#!env python3

import os, re, shutil
import eHive
import gzip
import csv, json

from Bio import SeqIO, SeqRecord
import requests, sys

class process_seq_region(eHive.BaseRunnable):

    def param_defaults(self):
        return {

                "synonym_map" : {
                    "Sequence-Name" : "INSDC_submitted_name",
                    "GenBank-Accn" : "INSDC",
                    "RefSeq-Accn" : "RefSeq",
                    "Assigned-Molecule" : "GenBank",
                    },

                "molecule_location" : {
                    "chromosome" : "nuclear_chromosome",
                    "mitochondrion" : "mitochondrial_chromosome",
                    "apicoplast" : "apicoplast_chromosome",
                    "plasmid" : "plasmid"
                    },
                "location_codon" : {
                    "apicoplast_chromosome" : 4
                    }
                }
        

    def run(self):
        genome_data = self.param('genome_data')
        work_dir = self.param('work_dir')
        report_path = self.param('report')
        gbff_path = self.param('gbff')

        # Create dedicated work dir
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final file name
        metadata_type = "seq_region"
        new_file_name = metadata_type + ".json"
        final_path = os.path.join(work_dir, new_file_name)

        # Get seq_regions data from report and gff3, and merge them
        report_regions = self.get_report_regions(report_path)
        gbff_regions = self.get_gbff_regions(gbff_path)
        seq_regions = self.merge_regions(report_regions, gbff_regions)
        
        # Setup the BRC4_seq_region_name
        self.add_brc4_ebi_name(seq_regions)

        # Guess translation table
        self.guess_translation_table(seq_regions)
        self.get_mitochondrial_codon_table(seq_regions, genome_data["species"]["taxonomy_id"])

        # Print out the file
        self.print_json(final_path, seq_regions)

        # Flow out the file and type
        output = {
                "metadata_type" : metadata_type,
                "metadata_json": final_path
                }
        self.dataflow(output, 2)

    def guess_translation_table(self, seq_regions) -> None:
        """
        Guess codon table based on location
        """
        location_codon = self.param("location_codon")
        
        for seqr in seq_regions:
            if "location" in seqr and seqr["location"] in location_codon:
                seqr["codon_table"] = location_codon[seqr["location"]]

    def add_brc4_ebi_name(self, seq_regions) -> None:
        """
        Use INSDC name without version as default BRC4 and EBI names
        """
        
        for seqr in seq_regions:
            if "synonyms" in seqr:
                for syn in seqr["synonyms"]:
                    if syn["source"] == "INSDC":
                        insdc_name = syn["name"]
                        flat_name, dot, version = insdc_name.partition(".")
                        seqr["BRC4_seq_region_name"] = flat_name
                        seqr["EBI_seq_region_name"] = flat_name

    def get_mitochondrial_codon_table(self, seq_regions, tax_id) -> None:
        """
        Get codon table for mitochondria based on taxonomy
        """
        tax_codon_table = self.get_mito_tax_codon_table(tax_id)
        if not tax_codon_table:
            return
        
        for seqr in seq_regions:
            if "location" in seqr and seqr["location"] == "mitochondrial_chromosome":
                if "codon_table" in seqr:
                    continue
                seqr["codon_table"] = tax_codon_table

    def get_mito_tax_codon_table(self, tax_id) -> int:
        """
        Get codon table for mitochondria based on taxonomy
        """
        server = "https://www.ebi.ac.uk"
        ext = "/ena/data/taxonomy/v1/taxon/tax-id/" + str(tax_id)

        r = requests.get(server + ext, headers={ "Content-Type" : "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        if "mitochondrialGeneticCode" in decoded:
            return int(decoded["mitochondrialGeneticCode"])
    
    def print_json(self, path, data) -> None:
        with open(path, "w") as json_out:
            json_out.write(json.dumps(data, sort_keys=True, indent=4))
    
    def merge_regions(self, regions1, regions2) -> list:
        """
        Merge seq_regions from different sources
        Return a list of seq_regions
        """
        if not regions1: regions1 = []
        if not regions2: regions2 = []
        
        # Get all the names
        names1 = frozenset(regions1)
        names2 = frozenset(regions2)
        all_names = names1.union(names2)
        
        # Create the seq_regions, merge if needed
        seq_regions = []
        for name in all_names:
            seqr1 = None
            seqr2 = None
            if name in names1:
                seqr1 = regions1[name]
            if name in names2:
                seqr2 = regions2[name]
            if seqr1 and seqr2:
                merged_seqr = {**seqr1, **seqr2}
                seq_regions.append(merged_seqr)
            elif seqr1:
                seq_regions.append(seqr1)
            else:
                seq_regions.append(seqr2)

        seq_regions.sort(key=lambda x: x["name"])
        return seq_regions
                

    def get_gbff_regions(self, gbff_path) -> dict:
        """
        Get seq_region data from the gbff file
        Return a dict of seq_regions, with their name as the key
        """
        if not gbff_path:
            return None
        
        seq_regions = {}
        _open = gbff_path.endswith(".gz") and gzip.open or open
        with _open(gbff_path, 'rt') as gbff_file:

            for record in SeqIO.parse(gbff_file, "genbank"):
                seqr = {}
                
                # Is the seq_region circular?
                annotations = record.annotations
                if "topology" in annotations and annotations["topology"] == "circular":
                    seqr["circular"] = True

                # Is there a genetic code defined?
                codon_table = self.get_codon_table(record)
                if codon_table:
                    seqr["codon_table"] = codon_table
                
                # Store the seq_region
                if seqr:
                    seq_regions[record.id] = seqr

        return seq_regions
    
    def get_codon_table(self, record) -> int:
        """
        Given a genbank record, seeks codon table features
        Returns a number if found
        """
        for feat in record.features:
            if "transl_table" in feat.qualifiers:
                return int(feat.qualifiers["transl_table"][0])
        
        return

    def get_report_regions(self, report_path) -> dict:
        """
        Get seq_region data from report file
        Return a dict of seq_regions, with their name as the key
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
            seq_region = self.make_seq_region(row, assembly_level)
            seq_regions[seq_region["name"]] = seq_region
        
        return seq_regions
    
    def make_seq_region(self, row, assembly_level) -> dict:
        """
        From a row of the report, create one seq_region
        Return a seq_region dict
        """
        seq_region = {}

        # Map the fields to their synonym name
        synonym_map = self.param("synonym_map")
        molecule_location = self.param("molecule_location")

        # Synonyms
        synonyms = []
        for field, source in synonym_map.items():
            if field in row and row[field].lower() != "na":
                synonym = { "source" : source, "name" : row[field] }
                synonyms.append(synonym)
        if len(synonyms) > 0:
            synonyms.sort(key=lambda x: x["source"])
            seq_region["synonyms"] = synonyms
        
        # Length
        field = "Sequence-Length"
        name = "length"
        if field in row and row[field].lower() != "na":
            seq_region[name] = int(row[field])
        
        # Name
        field = "GenBank-Accn"
        name = "name"
        if field in row and row[field].lower() != "na":
            seq_region[name] = row[field]
        
        # Coord system and location
        seq_role = row["Sequence-Role"]
        
        if seq_role in ("unplaced-scaffold", "unlocalized-scaffold"):
            seq_region["coord_system_level"] = assembly_level
        elif seq_role == "assembled-molecule":
            location = row["Assigned-Molecule-Location/Type"].lower()
            
            # Get location metadata
            if location in molecule_location:
                seq_region["coord_system_level"] = "chromosome"
                seq_region["location"] = molecule_location[location]
            else:
                raise Exception("Unrecognized sequence location: %s (is %s)" % (location, str(molecule_location)))
        else:
            raise Exception("Unrecognized sequence role: %s" % seq_role)
        return seq_region

    def report_to_csv(self, report_path) -> (str, dict):
        """
        Load an assembly report as a csv string
        Returns the csv as a string, and the head metadata as a dict
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
                        last_head = None
                    data += line
            return data, metadata

