#!env python3

import eHive
import json

class say_accession(eHive.BaseRunnable):

    def run(self):
        genome_data = self.param_required("genome_data")
        self.dataflow({ "accession" : genome_data["assembly"]["accession"] }, 2)
