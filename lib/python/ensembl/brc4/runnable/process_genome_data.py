#!env python3

import os, re, shutil
import eHive
import gzip
import csv, json

from Bio import SeqIO, SeqRecord

class process_genome_data(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                "provider" : {
                    "assembly" : {
                        "provider_name" : "genbank",
                        "provider_url" : "https://www.ncbi.nlm.nih.gov/assembly",
                        },
                    "annotation" : {
                        "provider_name" : "genbank",
                        "provider_url" : "https://www.ncbi.nlm.nih.gov/assembly",
                        },
                    }
                }
        

    def run(self):
        genome_data = self.param('genome_data')
        work_dir = self.param('work_dir')

        # Create dedicated work dir
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final file name
        metadata_type = "genome"
        new_file_name = metadata_type + ".json"
        final_path = os.path.join(work_dir, new_file_name)

        # Amend metadata
        self.add_provider(genome_data)

        # Print out the file
        self.print_json(final_path, genome_data)

        # Flow out the file and type
        output = {
                "metadata_type" : metadata_type,
                "metadata_json": final_path
                }
        self.dataflow(output, 2)
    
    def print_json(self, path, data) -> None:
        with open(path, "w") as json_out:
            json_out.write(json.dumps(data, sort_keys=True, indent=4))
    
    def add_provider(self, genome_data):
        """Add default provider metadata for assembly and gene models"""
        
        default_provider = self.param("provider")
        
        # RETROCOMPATIBILITY: move provider to assembly level
        if "provider" in genome_data:
            provider = genome_data["provider"]
            genome_data["assembly"]["provider_name"] = provider["name"]
            genome_data["assembly"]["provider_url"] = provider["url"]
            del genome_data["provider"]
        else:
            assembly = genome_data["assembly"]
            if not "provider_name" in assembly and not "provider_url" in assembly:
                assembly["provider_name"] = provider["assembly"]["provider_name"]
                assembly["provider_url"] = provider["assembly"]["provider_url"]
            genome_data["assembly"] = assembly
        
        # Annotation provider, if there are gene models
        if self.param("gff3_raw"):
            annotation = {}
            if "annotation" in genome_data:
                annotation = genome_data["annotation"]
            if not "provider_name" in annotation and not "provider_url" in annotation:
                annotation["provider_name"] = default_provider["annotation"]["provider_name"]
                annotation["provider_url"] = default_provider["annotation"]["provider_url"]
            genome_data["annotation"] = annotation
                
