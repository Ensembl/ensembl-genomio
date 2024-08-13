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


workflow ACCESSION_METADATA {
    // Generate a meta value with 2 keys: id and accession (same value)
    take:
        input_json

    emit:
        accession_meta
    
    main:
        accession_meta = _GET_ACCESSION_FROM_META_JSON(input_json).map {
            it, json -> tuple([id: it, accession: it], json)
        }
}

process _GET_ACCESSION_FROM_META_JSON {
    tag "$input_json"
    label 'local'

    input:
        path(input_json, stageAs: "input.json")

    output:
        tuple env(accession), path("input.json")
        
    shell:
        '''
        accession=$(jq -r '.["assembly"]["accession"]' !{input_json})
        '''
}
