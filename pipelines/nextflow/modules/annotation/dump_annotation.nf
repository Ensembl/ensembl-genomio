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

process DUMP_ANNOTATION {
    tag "${db.species}"
    label 'variable_2_8_32'
    maxForks params.max_database_forks

    input:
        val db

    output:
        tuple val(db), val("functional_annotation"), path("*.json")

    script:
        output = "${db.species}_functional_annotation.json"
        schema = "functional_annotation"
        password_arg = db.server.password ? "--pass ${db.server.password}" : ""
        """
        dump_functional_annotation.pl \
            --host $db.server.host \
            --port $db.server.port \
            --user $db.server.user \
            $password_arg \
            --dbname $db.server.database > $output
        schemas_json_validate --json_file $output --json_schema $schema
        """
    
    stub:
        output = "functional_annotation.json"
        """
        echo "[]" > $output
        """
}
