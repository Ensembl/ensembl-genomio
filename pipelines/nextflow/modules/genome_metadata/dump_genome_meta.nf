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


process DUMP_GENOME_META {
    tag "${db.species}"
    label 'normal'
    maxForks params.max_database_forks

    input:
        val db

    output:
        tuple val(db), val("genome"), path("*_genome.json")

    script:
        output = "${db.species}_genome.json"
        schema = "genome"
        password_arg = db.server.password ? "--password $db.server.password" : ""
        """
        genome_metadata_dump --host '$db.server.host' \
            --port '$db.server.port' \
            --user '$db.server.user' \
            $password_arg \
            --database '$db.server.database' \
            --verbose > $output
        schemas_json_validate --json_file $output --json_schema $schema
        """
    
    stub:
        output_file = "test_genome.json"
        schema = "genome"
        dump_dir = "$workflow.projectDir/../../../../data/test/pipelines/dumper/dump_files"
        dump_file = "dumped_genome.json"
        """
        cp $dump_dir/$dump_file $output_file
        schemas_json_validate --json_file $output_file --json_schema $schema
        """
}
