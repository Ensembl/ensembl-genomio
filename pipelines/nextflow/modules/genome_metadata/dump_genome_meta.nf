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
    label "normal"
    maxForks 10

    input:
        val server
        val db

    output:
        tuple val(db), val("genome"), path("genome.json")

    script:
        """
        genome_metadata_dump --host '${server.host}' \
            --port '${server.port}' \
            --user '${server.user}' \
            --password '${server.password}' \
            --database '${db.database}' \
            --debug \
            > genome.json
        """
}
