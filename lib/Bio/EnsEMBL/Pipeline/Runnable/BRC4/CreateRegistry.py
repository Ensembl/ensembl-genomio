#!env python3

import eHive
import re
import sys

from os import path

class CreateRegistry(eHive.BaseRunnable):

    def param_defaults(self):
        return {
        }

    def run(self):
        dbnames = self.param_required("db_name")
        registry_path = self.param_required("registry_path")
        
        if dbnames:
            print("Write to " + registry_path)
            
            with open(registry_path, "w") as registry:
                self.print_intro(registry)
                for dbname in dbnames:
                    self.print_db(registry, dbname)
                registry.write("\n1;\n")
    
    def print_intro(self, registry):
        
        intro = """use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

"""
        registry.write(intro)
        
    
    def print_db(self, registry, dbname):
        
        species = dbname.split("_core_", 1)[0]
        
        if self.param_exists('db_prefix'):
            db_prefix = self.param('db_prefix')
            if db_prefix:
                db_prefix = db_prefix + "_"
                species = species.split(db_prefix, 1)[1]
        
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

