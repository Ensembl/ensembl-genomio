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

    when:
        "seq_regions" in db.dump_selection

    script:
        output = "${db.species}_seq_region.json"
        schema = params.json_schemas["seq_region"]
        """
        seq_region_dump --host '$db.server.host' \
            --port '$db.server.port' \
            --user '$db.server.user' \
            --password '$db.server.password' \
            --database '$db.server.database' \
            > $output
        schemas_json_validate --json_file $output --json_schema $schema
        """
}
