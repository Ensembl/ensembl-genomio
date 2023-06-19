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

nextflow.enable.dsl=2
params.merge_split_genes="False"

process PROCESS_GFF3 {
    tag "$gff3 - $task.attempt"
    label 'variable_2_8_32'
    errorStrategy 'ignore'

    input:
        tuple val(gca), path(gff3), path(genome)

    output:
        tuple val(gca), path("*.json"), emit: functional_annotation
        tuple val(gca), path("*.gff3"), emit: gene_models

    script:
    """
    process_gff3 --in_gff_path ${gff3} --genome_data ${genome}
    """
}
