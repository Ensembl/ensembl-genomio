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

process DUMP_FASTA_PEPTIDES {
    tag "${db.species}"
    label "variable_2_8_32"
    maxForks params.max_database_forks

    input:
        val server
        val db
        val do_dump
    
    when:
        do_dump

    output:
        tuple val(db), val("fasta_pep"), path("*.fasta")

    script:
        output = "${db.species}_fasta_pep.fasta"
        """
        perl ${params.ensembl_root_dir}/ensembl-analysis/scripts/protein/dump_translations.pl \
            -host ${server.host} \
            -port ${server.port} \
            -user ${server.user} \
            -dbname ${db.database} \
            -dnadbhost ${server.host} \
            -dnadbport ${server.port} \
            -dnadbuser ${server.user} \
            -dnadbname ${db.database} \
            -file $output
        """
}
