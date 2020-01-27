#!env python3

import eHive
import json
from os import path

class JsonSchemaFactory(eHive.BaseRunnable):

    def run(self):
        manifest_path = self.param_required("manifest")

        errors = []

        with open(manifest_path) as manifest_file:
            manifest = json.load(manifest_file)

            # Check the manifest schema itself
            self.dataflow_json("manifest", manifest_path)
            
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

            # Check the other jsons schemas
            for metadata_type in ("seq_region", "genome", "functional_annotation"):
                if metadata_type in manifest:
                    self.dataflow_json(metadata_type, manifest[metadata_type])

    def dataflow_json(self, metadata_type, metadata_json):
        self.dataflow(
                {
                    "metadata_type": metadata_type,
                    "metadata_json": metadata_json
                    }, 2)
