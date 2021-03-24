#!env python3

import os, re, shutil
import eHive
import gzip
import csv, json
import datetime

from Bio import SeqIO, SeqRecord
import requests
import xml.etree.ElementTree as ET

class process_genome_data(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                "provider" : {
                    'GenBank' : {
                        "assembly" : {
                            "provider_name" : "GenBank",
                            "provider_url" : "https://www.ncbi.nlm.nih.gov/assembly",
                            },
                        "annotation" : {
                            "provider_name" : "GenBank",
                            "provider_url" : "https://www.ncbi.nlm.nih.gov/assembly",
                            },
                        },
                    'RefSeq' : {
                        "assembly" : {
                            "provider_name" : "RefSeq",
                            "provider_url" : "https://www.ncbi.nlm.nih.gov/refseq",
                            },
                        "annotation" : {
                            "provider_name" : "RefSeq",
                            "provider_url" : "https://www.ncbi.nlm.nih.gov/refseq",
                            },
                        },
                    },
                'accession_api_url' : "https://www.ebi.ac.uk/ena/browser/api/xml/%s",
                }
        

    def run(self):
        json_path = self.param_required("json_path")
        
        genome_data = self.get_json(json_path)

        # Amend metadata
        self.add_provider(genome_data)
        self.add_assembly_version(genome_data)
        self.add_genebuild_metadata(genome_data)
        self.add_species_metadata(genome_data)
        
        # Create dedicated work dir
        accession = genome_data["assembly"]["accession"]
        self.param("accession", accession)
        work_dir = self.param('work_dir')
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final file name
        metadata_type = "genome"
        new_file_name = metadata_type + ".json"
        final_path = os.path.join(work_dir, new_file_name)

        # Print out the file
        self.print_json(final_path, genome_data)

        # Flow out the file and type
        output = {
                "genome_json" : final_path,
                "genome_data" : genome_data,
                "accession" : accession
                }
        self.dataflow(output, 2)
    
    def get_json(self, json_path) -> dict:
        with open(json_path) as json_file:
            data = json.load(json_file)
            return data
    
    def print_json(self, path, data) -> None:
        with open(path, "w") as json_out:
            json_out.write(json.dumps(data, sort_keys=True, indent=4))
    
    def add_provider(self, genome_data):
        """Add default provider metadata for assembly and gene models"""
        
        # Provider = GenBank or RefSeq
        accession = genome_data["assembly"]["accession"]
        provider_data = self.param("provider")
        if accession.startswith("GCF"):
            provider = provider_data["RefSeq"]
        elif accession.startswith("GCA"):
            provider = provider_data["GenBank"]
        else:
            raise Exception("Accession doesn't look like an INSDC or RefSeq accession: " + accession)

        # Assembly provider
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
                annotation["provider_name"] = provider["annotation"]["provider_name"]
                annotation["provider_url"] = provider["annotation"]["provider_url"]
            genome_data["annotation"] = annotation
             
    def add_assembly_version(self, genome_data):
        """Add version number to assembly"""
        
        assembly = genome_data["assembly"]
        
        if not "version" in assembly:
            accession = assembly["accession"]
            values = accession.split(".")
            if len(values) == 2 and values[1]:
                assembly["version"] = int(values[1])

    def add_genebuild_metadata(self, genome_data):
        """Add metadata to genebuild"""
        
        assembly = genome_data["assembly"]
        genebuild = genome_data["genebuild"]
        
        current_date = datetime.date.today().isoformat()
        if not "version" in genebuild:
            genebuild["version"] = current_date
        if not "start_date" in genebuild:
            genebuild["start_date"] = current_date

    def add_species_metadata(self, genome_data):
        """Add species metadata from the accession"""
        
        species = genome_data["species"]
        accession = genome_data["assembly"]["accession"]
        
        if not "taxonomy_id" in species:
            taxonomy = self.get_taxonomy_from_accession(accession)
            species["taxonomy_id"] = taxonomy["taxon_id"]
            
            if (not "strain" in species) and "strain" in taxonomy:
                species["strain"] = taxonomy["strain"]
            
            if not"scientific_name" in species:
                species["scientific_name"] = taxonomy["scientific_name"]
    
    def get_taxonomy_from_accession(self, accession):
        """Provided an accession, get the associated taxonomy metadata"""
        
        url = self.param("accession_api_url")
        
        response = requests.get(url % accession)
        entry_xml = response.text
        
        entry = ET.fromstring(entry_xml)
        taxon_node = entry.find(".//TAXON")
        
        taxon_id = self.get_node_text(taxon_node, "TAXON_ID")
        strain = self.get_node_text(taxon_node, "STRAIN")
        scientific_name = self.get_node_text(taxon_node, "SCIENTIFIC_NAME")
        
        if not taxon_id:
            raise Exception("No taxon_id found for accession %s" % accession)
        if not scientific_name:
            raise Exception("No scientific_name found for accession %s" % accession)
        
        taxonomy = {
                'taxon_id' : int(taxon_id),
                'scientific_name' : scientific_name,
                }
        if strain:
            taxonomy['strain'] = strain
        
        return taxonomy
    
    def get_node_text(self, node, tag):
        try:
            return node.find(tag).text
        except:
            return
