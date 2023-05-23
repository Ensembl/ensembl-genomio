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
    label 'adaptive'
    tag "$gca - $task.attempt"
    debug true

    input:
    path genome_json
    path asm_report
    path genbank_gbff
    val gca
    val brc4_mode

    output:
    val gca, emit: gca
    path "${gca}/genome_amended.json", emit: amended_json
   
    script:
    """
    amend_genomic_data --genome_infile ${genome_json} --INSDC_RefSeq_report_infile ${asm_report} --genbank_infile ${genbank_gbff} --output_dir ${gca} --brc4_mode ${brc4_mode}
    """
}
