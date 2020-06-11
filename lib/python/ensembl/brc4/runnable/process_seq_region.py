#!env python3

import os
import eHive
import gzip
import shutil
import csv
import json

class process_seq_region(eHive.BaseRunnable):

    def run(self):
        genome_data = self.param('genome_data')
        work_dir = self.param('work_dir')
        report_path = self.param('report')
        gff3_path = self.param('gff3')

        # Set and create dedicated work dir
        accession = genome_data["assembly"]["accession"]
        work_dir = os.path.join(work_dir, accession)
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final file name
        metadata_type = "seq_region"
        new_file_name = accession + "_" + metadata_type + ".json"
        final_path = os.path.join(work_dir, new_file_name)

        # Extract the seq_region informations from the report
        report_regions = self.get_report_regions(report_path)

        # Get metadata from the gff3
        gff3_regions = self.get_gff3_regions(gff3_path)
        
        # Merge the metadata from gff3
        seq_regions = self.merge_regions(report_regions, gff3_regions)

        # Print out the file
        self.print_json(final_path, seq_regions)
        print(final_path)

        # No other operation
        output = {
                "metadata_type" : metadata_type,
                "metadata_json": final_path
                }
        self.dataflow(output, 2)
    
    def print_json(self, path, data) -> None:
        with open(path, "w") as json_out:
            json_out.write(json.dumps(data, sort_keys=True, indent=4))

    def report_to_csv(self, report_path) -> str:
        """Load an assembly report as a csv string"""

        _open = report_path.endswith(".gz") and gzip.open or open
        with _open(report_path, 'rt') as report:
            data = ""
            last_head = ""
            for line in report:
                # Ignore header
                if line.startswith("#"):
                    last_head = line
                    continue
                else:
                    if last_head:
                        data += last_head[2:].strip() + "\n"
                        last_head = None
                    data += line
            return data
    
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
                

    def get_gff3_regions(self, gff3_path) -> dict:
        """
        Get seq_region data from the gff3 if available
        Return a dict of seq_regions, with their name as the key
        """
        if not gff3_path:
            return None
        
        seq_regions = {}
        return seq_regions

    def get_report_regions(self, report_path) -> dict:
        """
        Get seq_region data from report file
        Return a dict of seq_regions, with their name as the key
        """

        # Map the fields to their synonym name
        synonym_map = {
                "Sequence-Name" : "INSDC_submitted_name",
                "GenBank-Accn" : "INSDC",
                "RefSeq-Accn" : "RefSeq",
                "Assigned-Molecule" : "GenBank",
                }

        # Get the report in a CSV format, easier to manipulate
        report_csv = self.report_to_csv(report_path)
        
        # Feed the csv string to the CSV reader
        reader = csv.DictReader(report_csv.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE)
        
        # Create the seq_regions
        seq_regions = {}
        for row in reader:
            seq_region = {}
            
            # Synonyms
            synonyms = []
            for field, source in synonym_map.items():
                if field in row and row[field].lower() != "na":
                    synonym = { "source" : source, "name" : row[field] }
                    synonyms.append(synonym)
            if len(synonyms) > 0:
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
            
            if seq_role == "unplaced-scaffold":
                seq_region["coord_system_level"] = "scaffold"
            elif seq_role == "assembled-molecule":
                location = row["Assigned-Molecule-Location/Type"]
                if location == "Mitochondrion":
                    seq_region["coord_system_level"] = "chromosome"
                    seq_region["location"] = "mitochondrial_chromosome"
                elif location == "Plasmid":
                    seq_region["coord_system_level"] = "chromosome"
                    seq_region["location"] = "plasmid"
                else:
                    raise Exception("Unrecognized sequence location: %s" % seq_location)
            else:
                raise Exception("Unrecognized sequence role: %s" % seq_role)
            
            seq_regions[seq_region["name"]] = seq_region
        
        return seq_regions
    
