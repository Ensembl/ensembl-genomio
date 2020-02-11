#!env python3

import eHive
import json
import mysql.connector
import re
import sys

from os import path


class PrepareGenome(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                "ensembl_mode" : False,
                "db_prefix" : ""
        }

    def run(self):
        manifest_path = self.param_required("manifest")

        errors = []
        
        manifest = self.get_manifest(manifest_path)
        genome = self.get_genome(manifest)
        self.update_assembly_name(genome)

        print(manifest)
        print("")
        print(genome)

        db_name = self.make_db_name(genome)
        print("db_name = %s" % db_name)
        self.param("db_name", db_name)

        # FLOW
        self.dataflow(
                {
                    "db_name" : db_name,
                    "manifest_data": manifest,
                    "genome_data": genome
                    }, 2)


    def get_manifest(self, manifest_path):
        with open(manifest_path) as manifest_file:
            manifest = json.load(manifest_file)
            
            # Use dir name from the manifest
            for name in manifest:
                if "file" in manifest[name]:
                    file_name = manifest[name]["file"]
                    manifest[name] = path.join(path.dirname(manifest_path), file_name)
                else:
                    for f in manifest[name]:
                        if "file" in manifest[name][f]:
                            file_name = manifest[name][f]["file"]
                            manifest[name][f] = path.join(path.dirname(manifest_path), file_name)
            return manifest
            
    def get_genome(self, manifest):
        # Get genome metadata
        if "genome" in manifest:
            genome_path = manifest["genome"]
            print(genome_path)

            with open(genome_path) as genome_file:
                return json.load(genome_file)

    def make_db_name(self, genome):
        
        production_name = self.make_production_name(genome)
        ensembl_version = str(self.param_required("ensembl_version"))
        release = str(self.param_required("release"))
        assembly_version = str(self.get_assembly_version(genome))
        db_prefix = self.param('db_prefix')

        name_parts = [
                production_name,
                "core",
                release,
                ensembl_version,
                assembly_version
                ]

        db_name = "_".join(name_parts)
        self.check_db_name_format(db_name)

        if db_prefix:
            db_name = db_prefix + "_" + db_name

        return db_name
    
    def make_production_name(self, genome):
        if "species" in genome:
            genome_sp = genome["species"]
            if "production_name" in genome_sp:
                prod_name = genome_sp["production_name"]
            else:
                # Create prod name from taxonomy
                if "taxonomy_id" in genome_sp:
                    taxon_id = genome_sp["taxonomy_id"]
                    scientific_name = self.get_scientific_name(taxon_id)
                    
                    # Only keep the first 2 words
                    split_name = scientific_name.split(" ")
                    genus = split_name[0].lower()
                    species = split_name[1]
                    
                    # Production name, with INSDC accession for uniqueness
                    accession = genome["assembly"]["accession"]
                    if not accession or accession == "":
                        raise Exception("The INSDC accession is needed")
                    accession = re.sub("\.\d+$", "", accession)
                    accession = accession.lower()
                    prod_name = genus + "_" + species + "_" + accession
                    genome["species"]["production_name"] = prod_name
                else:
                    raise Exception("Can't find taxonomy id for genome: %s" % genome)
            return prod_name
        else:
            raise Exception("Can't make production name for genome: %s" % genome)
    
    def get_scientific_name(self, taxon_id):
        
        host = self.param("taxonomy_host")
        user = self.param("taxonomy_user")
        port = self.param("taxonomy_port")
        dbname = self.param("taxonomy_dbname")
        
        # Get connector to taxon db
        mydb = mysql.connector.connect(
                host=host,
                user=user,
                port=port,
                database=dbname
                )
        cur = mydb.cursor()
        
        # Get scientific name
        select_st = "SELECT name FROM ncbi_taxa_name WHERE name_class=%s and taxon_id=%s"
        data = ("scientific name", taxon_id)
        cur.execute(select_st, data)
        row = cur.fetchone()
        for s in row:
            return s
        return

    def get_assembly_version(self, genome):

        if "assembly" in genome:
            assembly = genome["assembly"]
            if "version" in assembly:
                aversion = assembly["version"]
                if isinstance(aversion, int):
                    return aversion
                else:
                    try:
                        aversion = int(aversion)
                        genome["assembly"]["version"] = aversion
                        return aversion
                    except:
                        raise Exception("Assembly version is not an integer: %s" % aversion)
            elif "name" in assembly:
                # Get last number at the end of the assembly name
                matches = re.match("\w+?(\d+?)$", assembly["name"])
                if matches:
                    aversion = matches.group(1)
                    genome["assembly"]["version"] = aversion
                    return aversion
                else:
                    raise Exception("No assembly version found in %s" % str(genome))
            else:
                raise Exception("Can't get assembly version from name in from %s" % str(genome))
        else:
            raise Exception("Can't get assembly data in genome %s" % str(genome))
                

    def check_db_name_format(self, db_name):
        match = re.match("[a-z]+_[a-z]+(_[A-z0-9]+){0,3}_core_\d+_\d+_\d+$", db_name)

        if not match:
            raise Exception("Generated DB name is not the right format: %s" % db_name)


    def update_assembly_name(self, data):
        if data and "assembly" in data:
            asm = data["assembly"]
            if "name" not in asm:
                if "accession" not in asm:
                    raise Exception("no accession or name in genome_data/assembly")
                acc = asm["accession"].replace("_","").replace(".","v")
                asm["name"] = acc
                print("using \"%s\" as assembly.name" % (data["assembly"]["name"]), file = sys.stderr)
        else:
            raise Exception("no 'assembly' data provided in genome_data")
