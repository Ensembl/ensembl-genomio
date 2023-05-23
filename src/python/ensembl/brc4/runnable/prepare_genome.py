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


import json
from os import path
import re
import sys

import eHive
import mysql.connector


class prepare_genome(eHive.BaseRunnable):
    def param_defaults(self):
        return {"ensembl_mode": False, "db_prefix": ""}

    def run(self):
        manifest_path = self.param_required("manifest")

        errors = []

        manifest = self.get_manifest(manifest_path)
        genome = self.get_genome(manifest)
        self.update_accession(genome)
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
                "db_name": db_name,
                "manifest_data": manifest,
                "genome_data": genome,
                "species": genome["species"]["production_name"],
                "has_gff3": (1 if "gff3" in manifest and manifest["gff3"] else 0),
            },
            2,
        )

        # DB metadata for registry
        self.dataflow(
            {
                "db_data": {
                    "db_name": db_name,
                    "genome_data": genome,
                }
            },
            3,
        )

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
        prod_name = self.make_production_name(genome)
        prod_name = prod_name.replace(".", "")
        ensembl_version = str(self.param_required("ensembl_version"))
        release = str(self.param_required("release"))
        assembly_version = str(self.get_assembly_version(genome))
        db_prefix = self.param("db_prefix")

        # Check the db name is shorter than the limit: reduce it to "Genus sp"
        db_prod_name = prod_name
        max_prod_length = (
            64 - len("core") - len(str(release)) - len(str(ensembl_version)) - len(str(assembly_version)) - 4
        )
        if db_prefix:
            max_prod_length = max_prod_length - len(db_prefix) - 1
        if len(db_prod_name) > max_prod_length:
            print("DB name is too long! Trying to shorten it...")
            (genus_str, species_str, gca_str) = db_prod_name.split("_")
            db_prod_name = "_".join([genus_str, "sp", gca_str])

            # Still too long? Die
            if len(db_prod_name) > max_prod_length:
                raise Exception(
                    "Can't use the db name for %s (len = %d, max = %d)"
                    % (db_prod_name, len(db_prod_name), max_prod_length)
                )
            else:
                prod_name = db_prod_name

        # Store this production name
        genome["species"]["production_name"] = prod_name

        name_parts = [prod_name, "core", release, ensembl_version, assembly_version]

        db_name = "_".join(name_parts)
        db_name = db_name.replace(".", "")
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

                    # If what we get is not a species name (only 1 word), use the input name
                    if " " not in scientific_name:
                        print("Use input scientific name")
                        scientific_name = genome_sp["scientific_name"]

                    # Name in bracket means that it may be renamed by the community
                    # But we don't want those...
                    if "[" in scientific_name:
                        scientific_name = scientific_name.replace("[", "").replace("]", "")

                    if "-" in scientific_name:
                        scientific_name = scientific_name.replace("-", "")

                    # Only keep the first 2 words
                    split_name = scientific_name.split(" ")
                    genus = split_name[0].lower()
                    species = split_name[1].lower()

                    # Special case
                    if species == "sp.":
                        species = "sp"

                    # Check for special characters
                    if not re.match(r"^[a-zA-Z0-9]+$", genus):
                        raise Exception("The genus has special characters: " + genus)
                    if re.match(r"^[a-zA-Z0-9]$", species):
                        raise Exception("The species has special characters: " + species)
                    print(genus + " " + species)

                    # Production name, with INSDC accession for uniqueness
                    accession = genome["assembly"]["accession"]
                    if not accession or accession == "":
                        raise Exception("The INSDC accession is needed")
                    accession = re.sub("\.\d+$", "", accession).replace("_", "")
                    accession = accession.lower()
                    prod_name = genus + "_" + species + "_" + accession
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
        mydb = mysql.connector.connect(host=host, user=user, port=port, database=dbname)
        cur = mydb.cursor()

        # Get scientific name
        select_st = "SELECT name FROM ncbi_taxa_name WHERE name_class=%s and taxon_id=%s"
        data = ("scientific name", taxon_id)
        cur.execute(select_st, data)
        row = cur.fetchone()

        if row:
            for s in row:
                return s
        else:
            return

    def get_assembly_version(self, genome):
        if "assembly" in genome:
            assembly = genome["assembly"]
            if "version" in assembly:
                aversion = assembly["version"]

                if str(aversion) == "0":
                    raise Exception("Assembly version cannot be 0")

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
        match = re.match("[a-z]+_[a-z]+(_[A-z0-9]+){0,2}_core_\d+_\d+_\d+$", db_name)

        if not match:
            raise Exception(
                "Generated DB name is not the right format: %s (in %s)"
                % (db_name, self.param_required("manifest"))
            )

    def update_assembly_name(self, data):
        if data and "assembly" in data:
            asm = data["assembly"]
            if "name" not in asm:
                if "accession" not in asm:
                    raise Exception("no accession or name in genome_data/assembly")
                acc = asm["accession"].replace("_", "").replace(".", "v")
                asm["name"] = acc
                print('using "%s" as assembly.name' % (data["assembly"]["name"]), file=sys.stderr)
        else:
            raise Exception("no 'assembly' data provided in genome_data")

    def update_accession(self, data):
        if data and "assembly" in data:
            asm = data["assembly"]
            if "accession" in asm and asm["accession"].startswith("GCF_"):
                ### CHANGE TO GCA
                data["assembly"]["accession"] = data["assembly"]["accession"].replace("GCF_", "GCA_")
                print('using "%s" as assembly.accession' % (data["assembly"]["accession"]), file=sys.stderr)
        else:
            raise Exception("no 'assembly' data provided in genome_data")
