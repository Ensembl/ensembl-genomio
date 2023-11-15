
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

process PROCESS_FASTA {
    tag "${meta.accession}"
    label 'adaptive'

    input:
        tuple val(meta), path(compressed_gff), path(fasta_file), path(gbff_file)
        val pep_mode

    output:
        tuple val(meta), path ("fasta_*.fa"), emit: processed_fasta

    shell:
    output_fasta = pep_mode == 1 ? "fasta_pep.fa" : "fasta_dna.fa"
    pep_mode_flag = pep_mode == 0 ? "" : "--peptide_mode"
    '''
    fasta_process --fasta_infile !{fasta_file} \
        --genbank_infile !{gbff_file} \
        --fasta_outfile !{output_fasta} \
        !{pep_mode_flag}
    '''
}
