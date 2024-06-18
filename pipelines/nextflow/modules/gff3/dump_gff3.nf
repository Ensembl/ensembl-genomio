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

process DUMP_GFF3 {
    tag "${db.species}"
    label "variable_2_8_32"
    label "ensembl_scripts_container"
    maxForks params.max_database_forks

    input:
        val db

    output:
        tuple val(db), val("gff3"), path("*.gff3")

    script:
        output = "${db.species}.gff3"
        """
        dump_gff3.pl \
            --host $db.server.host \
            --port $db.server.port \
            --user $db.server.user \
            --pass $db.server.password \
            --dbname $db.server.database > $output
        """
    
    stub:
        output = "gene_models.gff3"
        """
        echo "##gff-version 3\n" > $output
        """
}
