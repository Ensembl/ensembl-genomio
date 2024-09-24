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


process PROCESS_GFF3 {
    tag "${meta.accession}"
    label 'adaptive'

    input:
        tuple val(meta), path(genome), path(gff3)

    output:
        tuple val(meta), path("functional_annotation.json"), emit: functional_annotation
        tuple val(meta), path("gene_models.gff3"), emit: gene_models

    shell:
        out_func = "functional_annotation.json"
        out_gff = "gene_models.gff3"
        '''
        gff3_process --genome_data !{genome} --in_gff_path !{gff3} --out_gff_path !{out_gff} \
            --out_func_path !{out_func} -v \
            --log_file gff3_process.log
        
        schemas_json_validate --json_file !{out_func} --json_schema functional_annotation
        '''
}
