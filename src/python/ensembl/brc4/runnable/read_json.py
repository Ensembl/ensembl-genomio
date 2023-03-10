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


import eHive
import json

class read_json(eHive.BaseRunnable):
    """Read a json data from a file, and flow it out for the pipeline to use.
    
    Args:
        json_path: json to load
        name: key of the dataflow object associated with the json data
    
    Dataflows:
        2: one Dict record with the name as key, and the json data as value
    """

    def run(self):
        json_path = self.param_required("json_path")
        name = self.param_required("name")
        
        with open(json_path) as json_file:
            data = json.load(json_file)
            self.dataflow({ name : data }, 2)
