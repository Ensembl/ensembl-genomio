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


process DUMP_NCBI_STATS {
    tag "${db.species}"
    label 'local'
    local 'cached'

    input:
        val server
        val db

    output:
        tuple val(db), path("ncbi_stats.json")

    shell:
        output = "ncbi_stats.json"
        '''
        function get_meta_value {
            meta_key=$1

            mysql --host="!{server.host}" --port="!{server.port}" --user="!{server.user}" \
                --password="!{server.password}" --database="!{db.database}" \
                -N -e "SELECT meta_value FROM meta WHERE meta_key='$meta_key'"
        }

        touch !{output}
        # Get the INSDC accession to use
        accession=$(get_meta_value "assembly.accession")
        provider=$(get_meta_value "assembly.provider_url" | sed -r 's%^.+/%%g')
        if [ $provider == "refseq" ]; then
            accession=$(echo $accession | sed 's/^GCA_/GCF_/')
        fi
        echo "Provider is $provider"
        echo "Accession is $accession"

        datasets summary genome accession $accession | jq '.' > !{output}

        # Check if it should maybe be using RefSeq?           
        if [ "$(jq '.total_count' !{output})" == "0" ]; then
            accession=$(echo $accession | sed 's/^GCA_/GCF_/')
            echo "Trying again with accession $accession"
            datasets summary genome accession $accession | jq '.' > !{output}
        fi
        '''
}