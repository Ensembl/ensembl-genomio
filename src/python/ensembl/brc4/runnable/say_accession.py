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

class say_accession(eHive.BaseRunnable):
    """Simple runnable to bring out the accession value for the pipeline to use.
    
    Args:
        genome_data: a dict from genome_data following the schema from schemas/genome_schema.json
    
    Dataflows:
        2: a single value named accession
    """

    def run(self):
        genome_data = self.param_required("genome_data")
        accession = genome_data["assembly"]["accession"]
        self.dataflow({ "accession" : accession }, 2)
