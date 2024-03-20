// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

process DATASETS_METADATA {
    label 'local'
    label 'cached'

    input:
        path(input_json, stageAs: "input_genome.json")

    output:
        tuple path("input_genome.json", includeInputs: true), path("ncbi_meta.json")
        
    shell:
        '''
        accession=$(jq -r '.["assembly"]["accession"]' !{input_json})
        datasets summary genome accession $accession | jq '.' > ncbi_meta.json
        if [ ! -s "ncbi_meta.json" ]; then
            echo "No Metadata from datasets for $accession"
            exit 1
        fi
        '''
    
    stub:
        """
        echo '{"reports":[{"organism":{"tax_id":1000}}]}' > ncbi_meta.json
        """
}
