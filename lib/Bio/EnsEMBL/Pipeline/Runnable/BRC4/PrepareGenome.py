#!env python3

import eHive
import json
from os import path
import re

class PrepareGenome(eHive.BaseRunnable):

    def param_default(self):
        return {
                ensembl_mode : False,
                db_prefix: ""
        }

    def run(self):
        manifest_path = self.param_required("manifest")

        errors = []
        
        manifest = self.get_manifest(manifest_path)
        genome = self.get_genome(manifest)

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

        production_name = genome["species"]["production_name"]
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

    def get_assembly_version(self, genome):

        if "assembly" in genome:
            assembly = genome["assembly"]
            if "version" in assembly and isinstance(assembly["version"], int):
                return assembly["version"]
            elif "name" in assembly:
                # Get last number at the end of the assembly name
                matches = re.match("\w+?(\d+?)$", assembly["name"])
                if matches:
                    version = matches.group(1)
                    return version

    def check_db_name_format(self, db_name):
        match = re.match("[a-z]+_[a-z]+(_[a-z0-9]+)?_core_\d+_\d+_\d+$", db_name)

        if not match:
            raise Exception("Generated DB name is not the right format: %s" % db_name)

