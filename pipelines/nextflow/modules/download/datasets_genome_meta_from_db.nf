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


process DOWNLOAD_GENOME_META_FROM_DB {
    tag "${db.species}"
    label 'local'
    label 'cached'
    label 'datasets_container'

    input:
        tuple val(db), val(accession)

    output:
        tuple val(db), path("ncbi_stats.json")

    shell:
        output = "ncbi_stats.json"
        password_arg = db.server.password ? "--password $db.server.password" : ""
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
        output_file = "ncbi_stats.json"
        dump_dir = "$workflow.projectDir/../../../../data/test/pipelines/dumper/dump_files"
        dump_file = "downloaded_ncbi_stats.json"
        """
        cp $dump_dir/$dump_file $output_file
        """
}
