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

process DUMP_SEQ_REGIONS {
    tag "${db.species}"
    label "variable_2_8_32"
    maxForks params.max_database_forks

    input:
        val db
    
    output:
        tuple val(db), val("seq_region"), path("*_seq_region.json")

    script:
        output = "${db.species}_seq_region.json"
        schema = "seq_region"
        password_arg = db.server.password ? "--password $db.server.password" : ""
        """
        seq_region_dump --host '$db.server.host' \
            --port '$db.server.port' \
            --user '$db.server.user' \
            $password_arg \
            --database '$db.server.database' \
            > $output
        schemas_json_validate --json_file $output --json_schema $schema
        """
    
    stub:
        output_file = "test_seq_region.json"
        schema = "seq_region"
        dump_dir = "$workflow.projectDir/../../../../data/test/pipelines/dumper/dump_files"
        dump_file = "dumped_seq_region.json"
        """
        cp $dump_dir/$dump_file $output_file
        schemas_json_validate --json_file $output_file --json_schema $schema
        """
}
