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


process DOWNLOAD_NCBI_STATS {
    tag "$meta.id"
    label 'local'
    label 'cached'
    container 'ensemblorg/datasets-cli:latest'

    input:
        val(meta)  // with keys [ id, accession ]

    output:
        tuple val(meta), path("ncbi_stats.json")

    shell:
        output = "ncbi_stats.json"
        '''
        echo "Calling datasets-cli.... datasets 'summary' 'genome' 'accession' [!{meta.accession}]'"

        # Pipe datasets to jq instead of '--as-json-lines' to 
        # obtain a total_count of reports returned.
        datasets summary genome accession !{meta.accession} | jq '.' > !{output}

        if [ "$?" -ne 0 ]; then
            echo "Invalid or unsupported assembly accession: !{meta.accession}"
            exit 1
        fi

        # Check if it should maybe be using RefSeq?
        if [[ $(jq '.total_count' !{output}) -eq 0 ]] && [[ !{meta.accession} =~ "GCA_" ]]; then
            accession=$(echo !{meta.accession} | sed 's/^GCA_/GCF_/')
            echo "Trying again with RefSeq accession: $accession"
            datasets summary genome accession $accession | jq '.' > !{output}
        fi
        '''
    
    stub:
        output_file = "ncbi_stats.json"
        dump_dir = "$workflow.projectDir/../../../../data/test/pipelines/dumper/dump_files"
        dump_file = "downloaded_ncbi_stats.json"
        """
        cp $dump_dir/$dump_file $output_file
        """
}
