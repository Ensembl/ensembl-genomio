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


process CORE_TO_ASM_META {
    tag "${db.species}"
    label 'local'
    label 'cached'

    input:
        val db

    output:
        tuple val(db), path("accession_meta.json")

    shell:
        output = "accession_meta.json"
        password_arg = db.server.password ? "--password $db.server.password" : ""
        '''
        function get_meta_value {
            meta_key=$1

            mysql \
                --host="!{db.server.host}" \
                --port="!{db.server.port}" \
                --user="!{db.server.user}" \
                !{password_arg} \
                --database="!{db.server.database}" \
                -N -e "SELECT meta_value FROM meta WHERE meta_key='$meta_key'"
        }

        touch !{output}
        # Get the INSDC accession to use
        accession=$(get_meta_value "assembly.accession")
        provider=$(get_meta_value "assembly.provider_url" | sed -r 's%^.+/%%g')
        if [ $provider == "refseq" ]; then
            accession=$(echo $accession | sed 's/^GCA_/GCF_/')
        fi

        echo -e -n "Provider is $provider\nAccession is $accession"
        echo "{@Prodvider@: @$provider@,@Accession@: @$accession@}" | sed 's/@/"/g' > !{output}        
        '''
    
    stub:
        output_file = "accession_meta.tsv"
        dump_dir = "$workflow.projectDir/../../../../data/test/pipelines/dumper/dump_files"
        dump_file = "dumped_core_accn_meta.json"
        """
        cp $dump_dir/$dump_file $output_file
        """
}
