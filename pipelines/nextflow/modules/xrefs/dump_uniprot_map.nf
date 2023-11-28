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

process DUMP_UNIPROT_MAP {
    tag "${db.species}"
    label "normal"
    maxForks params.max_database_forks

    input:
        val db
    
    output:
        path("*_uniprot_map.tsv")

    script:
        output = "${db.species}_uniprot_map.tsv"
        """
        touch $output
        xrefs_dump_uniprot_map \
            --host '${db.server.host}' \
            --port '${db.server.port}' \
            --user '${db.server.user}' \
            --password '${db.server.password}' \
            --database '${db.server.database}' \
            --verbose \
            > "$output"
        exit
        """
    
    stub:
        output = "${db.species}_uniprot_map.tsv"
        """
        touch $output
        """
}
