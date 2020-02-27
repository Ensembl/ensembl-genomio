#!env python3

import eHive
import json
from jsonschema import validate

class SchemaValidator(eHive.BaseRunnable):

    def param_defaults(self):
        return {
        }

    def run(self):
        json_file = self.param_required("json_file")
        json_schemas = self.param_required("json_schema")
        metadata_type = self.param_required("metadata_type")
        
        if metadata_type in json_schemas:
            json_schema = json_schemas[metadata_type]
        else:
            raise Exception("Schema not defined: %s" % metadata_type)
        
        file = self.get_json(json_file)
        schema = self.get_json(json_schema)
        
        validate(instance=file, schema=schema)

    def get_json(self, json_path):
        with open(json_path) as json_file:
            data = json.load(json_file)
            return data

