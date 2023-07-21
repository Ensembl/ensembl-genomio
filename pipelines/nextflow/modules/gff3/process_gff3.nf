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
    tag "$gca"
    label 'adaptive'

    input:
        tuple val(gca), path(gff3), path(genome)

    output:
        tuple val(gca), path("functional_annotation.json"), emit: functional_annotation
        tuple val(gca), path("gene_models.gff3"), emit: gene_models

    script:
    def out_func = "functional_annotation.json"
    def out_gff = "gene_models.gff3"
    """
    process_gff3 --genome_data $genome --in_gff_path $gff3 --out_gff_path $out_gff --out_func_path $out_func
    """
}
