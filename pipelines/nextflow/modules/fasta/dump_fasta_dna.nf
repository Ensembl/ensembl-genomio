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

process DUMP_FASTA_DNA {
    tag "${db.species}"
    label "variable_2_8_32"
    maxForks params.max_database_forks

    input:
        val db

    output:
        tuple val(db), val("fasta_dna"), path("*.fasta")

    script:
        output = "${db.species}_fasta_dna.fasta"
        """
        perl ${params.ensembl_root_dir}/ensembl-analysis/scripts/sequence_dump.pl \
            -dbhost ${db.server.host} \
            -dbport ${db.server.port} \
            -dbuser ${db.server.user} \
            -dbname ${db.server.database} \
            -coord_system_name toplevel \
            -toplevel \
            -onefile \
            -nonref \
            -filename $output
        """
    
    stub:
        output_file = "dna.fasta"
        files_dir = "$workflow.projectDir/../../../../data/test/pipelines/dumper/dump_files"
        """
        ln -s $files_dir/dumped_dna.fasta $output_file
        """
}
