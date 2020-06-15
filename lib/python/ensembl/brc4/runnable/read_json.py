#!env python3

import eHive
import json

class read_json(eHive.BaseRunnable):

    def run(self):
        json_path = self.param_required("json_path")
        name = self.param_required("name")
        
        with open(json_path) as json_file:
            data = json.load(json_file)
            self.dataflow({ name : data }, 2)
