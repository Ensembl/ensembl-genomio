#!env python3

import eHive
import re
import sys

from os import path

class CreateRegistry(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                "use_alias_name": False,
        }

    def run(self):
        db_data = self.param_required("db_data")
        registry_path = self.param_required("registry_path")
        
        if db_data:
            print("Write to " + registry_path)
            
            with open(registry_path, "w") as registry:
                self.print_intro(registry)
                for db in db_data:
                    dbname = db["db_name"]
                    prod_name = db["genome_data"]["species"]["production_name"]
                    
                    # Alias?
                    alias_name = prod_name
                    if "alias" in db["genome_data"]["species"]:
                        alias = db["genome_data"]["species"]["alias"]
                        if isinstance(alias, list):
                            alias_name = alias[0]
                        else:
                            alias_name = alias
                    
                    # Use prod_name or alias as name?
                    species = prod_name
                    if self.param("use_alias_name"):
                        species = alias_name
                    
                    self.print_db(registry, dbname, species)
                registry.write("\n1;\n")
    
    def print_intro(self, registry):
        
        intro = """use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

"""
        registry.write(intro)
        
    
    def print_db(self, registry, dbname, species):
        
        host = self.param("dbsrv_host")
        port = self.param("dbsrv_port")
        user = self.param("dbsrv_user")
        password = self.param("dbsrv_pass")
        
        adaptor = (
                "Bio::EnsEMBL::DBSQL::DBAdaptor->new(",
                "   '-species' => '%s'," % species,
                "   '-dbname' => '%s'," % dbname,
                "   '-group' => 'core',",
                "   '-host' => '%s'," % host,
                "   '-port' => '%d'," % port,
                "   '-user' => '%s'," % user,
                "   '-pass' => '%s'," % password,
                ");",
                ""
        )
        
        registry.write("\n".join(adaptor))

