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


from pathlib import Path

import eHive
from jsonschema import validate

from ensembl.brc4.runnable.utils import get_json


class schema_validator(eHive.BaseRunnable):
    """Check a json file with a provided json schema.

    Args:
        json_file: Path to the json to check.
        json_schema: Dict of json schema paths.
        metadata_type: Key to find to schema path to use from the json_schema.

    Dataflows:
        None
    """

    def run(self):
        json_file = Path(self.param_required("json_file"))
        json_schemas = self.param_required("json_schema")
        metadata_type = self.param_required("metadata_type")
        
        if metadata_type in json_schemas:
            json_schema = Path(json_schemas[metadata_type])
        else:
            raise Exception(f"Schema not defined: {metadata_type}")
        
        file = get_json(json_file)
        schema = get_json(json_schema)
        
        validate(instance=file, schema=schema)
