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

process AMEND_GENOME_DATA {
    tag "${meta.accession}"
    label 'adaptive'

    input:
        tuple val(meta), path(genome_json, stageAs: "incoming_genome.json"), path(asm_report),
            path(genomic_fna), path(genbank_gbff)

    output:
        tuple val(meta), path ("genome.json"), emit: amended_json
    
    shell:
        output = "genome.json"
        '''
        genome_metadata_extend --genome_infile !{genome_json} \
            --report_file !{asm_report} \
            --genbank_file !{genbank_gbff} \
            --genome_outfile !{output}
        
        schemas_json_validate --json_file !{output} --json_schema genome
        '''
}
