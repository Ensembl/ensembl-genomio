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

process DUMP_DB {
    publishDir "$out_dir/build_$db.release/coredb/$db.division", mode: 'copy'
    tag "$db.species"
    label "variable_2_8_32"
    maxForks params.max_database_forks

    input:
        val db
        val out_dir

    output:
        path "*.sql.gz"

    script:
        output_file = "${db.species}.sql.gz"
        """
        db_pass=""
        if [ "${db.server.password}" != "" ]; then
            db_pass="--password '${db.server.password}'"
        fi

        mysqldump '${db.server.database}' \
            --host '${db.server.host}' \
            --port '${db.server.port}' \
            --user '${db.server.user}' \
            \$db_pass \
            | gzip > $output_file
        """
    
    stub:
        output_file = "${db.species}.sql.gz"
        dump_dir = "$workflow.projectDir/../../../../data/test/pipelines/dumper/dump_sql"
        dump_file = "dumped_db.sql.test"
        """
        cat $dump_dir/$dump_file | gzip > $output_file
        """
}
