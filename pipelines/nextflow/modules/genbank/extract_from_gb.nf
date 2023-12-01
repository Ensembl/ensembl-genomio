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

process EXTRACT_FROM_GB {
    tag "$gb_file"
    label 'default'

    input:
        tuple val(meta), path(gb_file)

    output:
        tuple val(meta), path("genome.json"), emit: genome
        tuple val(meta), path("seq_region.json"), emit: seq_regions
        tuple val(meta), path("dna.fasta"), emit: dna_fasta
        tuple val(meta), path("*.gff"), emit: gene_gff, optional: true
        tuple val(meta), path("pep.fasta"), emit: pep_fasta, optional: true

    shell:
    '''
    genbank_extract_data \
        --prefix !{meta.prefix} \
        --prod_name !{meta.production_name} \
        --gb_file !{gb_file} \
        --debug

    schemas_json_validate --json_file "genome.json" --json_schema "genome"
    schemas_json_validate --json_file "seq_region.json" --json_schema "seq_region"
    '''
}
