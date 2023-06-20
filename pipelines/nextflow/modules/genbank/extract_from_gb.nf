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

process EXTRACT_FROM_GB {
    tag "$gb_file"
    label 'default'

    input:
        tuple val(accession), path(gb_file)
        val prefix
        val production_name

    output:
        tuple val(accession), path("*.gff"), emit: gene_gff
        tuple val(accession), path("genome.json"), emit: genome 
        tuple val(accession), path("seq_region.json"), emit: seq_regions
        tuple val(accession), path("dna.fasta"), emit: dna_fasta
        tuple val(accession), path("pep.fasta"), emit: pep_fasta, optional: true

    script:
    """
    extract_from_gb --prefix ${prefix} --prod_name ${production_name} --gb_file ${gb_file}
    """
}

