
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
    tag "${gca}"
    label 'adaptative'

    input:
    path fasta_file
    path gbff_file
    val gca
    val pep_mode

    output:
    val gca, emit: gca
    path "${gca}/*.fa", emit: processed_fasta

    script:
    """
    prep_fasta_data --fasta_infile ${fasta_file} --genbank_infile ${gbff_file} --output_dir "${gca}" --peptide_mode ${pep_mode}
    """
}
