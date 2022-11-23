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
import datetime
from pathlib import Path

import eHive
import requests
import xml.etree.ElementTree as ET

from ensembl.brc4.runnable.utils import Utils

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
        json_path = Path(self.param_required("json_path"))
        
        genome_data = Utils.get_json(json_path)

        # Amend metadata
        self.add_provider(genome_data)
        self.add_assembly_version(genome_data)
        self.add_genebuild_metadata(genome_data)
        self.add_species_metadata(genome_data)
        
        # Create dedicated work dir
        accession = genome_data["assembly"]["accession"]
        self.param("accession", accession)
        work_dir = Path(self.param('work_dir'))
        if not work_dir.is_dir():
            os.makedirs(work_dir)

        # Final file name
        metadata_type = "genome"
        new_file_name = metadata_type + ".json"
        final_path = work_dir / new_file_name

        # Print out the file
        Utils.print_json(final_path, genome_data)

        # Flow out the file and type
        output = {
                "genome_json" : str(final_path),
                "genome_data" : genome_data,
                "accession" : accession
                }
        self.dataflow(output, 2)
    
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
        
        # Use the GenBank accession without version
        gb_accession = accession.replace("GCF", "GCA")
        gb_accession = gb_accession.split(".")[0]

        response = requests.get(url % gb_accession)
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
