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

process DOWNLOAD_ASM_DATA {
    label 'adaptative'
    tag "$gca - $task.attempt"
    debug true

    input:
    val(gca)
    
    output:
    val gca, emit: gca
    path "${gca}/*_assembly_report.txt", emit: asm_report
    path "${gca}/*_genomic.fna.gz", emit: genome_fna
    path "${gca}/*_genomic.gbff.gz", emit: genome_gbff, optional: true
    path "${gca}/*_genomic.gff.gz", emit: genome_gff, optional: true
    path "${gca}/*_protein.faa.gz", emit: protein_fa, optional: true

    script:
    """
    retrieve_assembly_data --accession $gca --asm_download_dir ./
    """
}
