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

process DOWNLOAD_GENOME_META_FROM_ACC {
    tag "$accession"
    label 'local'
    label 'cached'
    label 'datasets_container'

    input:
        val(accession)

    output:
        tuple val(accession), path("ncbi_meta.json")
        
    shell:
        output = "ncbi_meta.json"
        '''
        echo "Calling datasets-cli.... datasets 'summary' 'genome' 'accession' [!{accession}]'"

        # Pipe datasets to jq instead of '--as-json-lines' to 
        # obtain a total_count of reports returned.
        datasets summary genome accession !{accession} | jq '.' > !{output}

        if [ "$?" -ne 0 ]; then
            echo "Invalid or unsupported assembly accession: !{accession}"
            exit 1
        fi

        # Check if it should maybe be using RefSeq?           
        if [[ $(jq '.total_count' !{output}) -eq 0 ]] && [[ !{accession} =~ "GCA_" ]]; then
            accession=$(echo !{accession} | sed 's/^GCA_/GCF_/')
            echo "Trying again with RefSeq accession: $accession"
            datasets summary genome accession !{accession} | jq '.' > !{output}
        fi
        '''
    
    stub:
        """
        echo '{"reports":[{"organism":{"tax_id":1000}}]}' > ncbi_meta.json
        """
}
