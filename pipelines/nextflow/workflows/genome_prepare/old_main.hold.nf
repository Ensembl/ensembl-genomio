#!/usr/bin/env nextflow
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

// Genome prepare pipeline
// See full documentation in docs/genome_prepare.md

// Print usage
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the 'Genome prepare' pipeline is as follows:

        nextflow run workflows/genome_prepare/main.nf --input_dir <json_input_dir>

        Mandatory arguments:
        --input_dir                    Location of input json(s) with component/organism genome metadata

        Optional arguments:
        --help                         This usage statement.
        """
}

// Check mandatory parameters
if (params.help) {
    helpMessage()
    exit 0
}
else if (params.input_dir) {
    ch_genome_json = Channel.fromPath("${params.input_dir}/*.json", checkIfExists: true)
} else {
    exit 1, 'Input directory not specified!'
}

// if (params.regions_to_exclude) {
//     params.regions_to_exclude_list = params.regions_to_exclude?.split(',') as List
// }


// Import modules/subworkflows
include { PREPARE_GENOME_METADATA } from '../../modules/prepare_genome_metadata.nf'
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_GENOME } from '../../modules/check_json_schema.nf'
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_SEQREG } from '../../modules/check_json_schema.nf'
include { DOWNLOAD_ASM_DATA } from '../../modules/download_asm_data.nf'
include { PROCESS_GFF3; GFF3_VALIDATION } from '../../modules/process_gff3.nf'
include { UNPACK_FILE } from '../../modules/unpack_gff3.nf'
include { PROCESS_SEQ_REGION } from '../../modules/process_seq_region.nf'
include { PROCESS_FASTA as PROCESS_FASTA_DNA } from '../../modules/process_fasta_data.nf'
include { PROCESS_FASTA as PROCESS_FASTA_PEP } from '../../modules/process_fasta_data.nf'
include { AMEND_GENOME_DATA } from '../../modules/amend_genome_data.nf'

// Run main workflow
workflow {
    ch_metadata_json = PREPARE_GENOME_METADATA(ch_genome_json)
    CHECK_JSON_SCHEMA_GENOME('genome', ch_metadata_json)
    DOWNLOAD_ASM_DATA(CHECK_JSON_SCHEMA_GENOME.out.gca)
    if(DOWNLOAD_ASM_DATA.out.gene_gff && DOWNLOAD_ASM_DATA.out.protein_fa){
        // println "GFF and Pep Present:"
        unpacked_gff = UNPACK_FILE(DOWNLOAD_ASM_DATA.out.gene_gff, 'gff')
        unpacked_gff.view()
        PROCESS_GFF3(unpacked_gff, CHECK_JSON_SCHEMA_GENOME.out.json_file)
        PROCESS_GFF3.out.functional_annotation.view()
        PROCESS_FASTA_PEP(DOWNLOAD_ASM_DATA.out.protein_fa, DOWNLOAD_ASM_DATA.out.genome_gbff, CHECK_JSON_SCHEMA_GENOME.out.gca, '1')
    }
    CHECK_JSON_SCHEMA_GENOME.out.gca.view()
    // DOWNLOAD_ASM_DATA.out.asm_report.view()
    // DOWNLOAD_ASM_DATA.out.genome_fna.view()
    // DOWNLOAD_ASM_DATA.out.asm_report.view()
    ch_process_seqregion = PROCESS_SEQ_REGION(CHECK_JSON_SCHEMA_GENOME.out.json_file, DOWNLOAD_ASM_DATA.out.asm_report, DOWNLOAD_ASM_DATA.out.genome_gbff, DOWNLOAD_ASM_DATA.out.gca)
    // ch_process_seqregion.view()
    CHECK_JSON_SCHEMA_SEQREG('seq_region', ch_process_seqregion)
    // CHECK_JSON_SCHEMA_SEQREG.out.gca.view()
    // CHECK_JSON_SCHEMA_SEQREG.out.json_file.view()
    PROCESS_FASTA_DNA(DOWNLOAD_ASM_DATA.out.genome_fna, DOWNLOAD_ASM_DATA.out.genome_gbff, CHECK_JSON_SCHEMA_SEQREG.out.gca, '0')
    PROCESS_FASTA_DNA.out.processed_fasta.view()
    AMEND_GENOME_DATA(CHECK_JSON_SCHEMA_GENOME.out.json_file, DOWNLOAD_ASM_DATA.out.asm_report, DOWNLOAD_ASM_DATA.out.genome_gbff, DOWNLOAD_ASM_DATA.out.gca, params.brc4_mode)
    AMEND_GENOME_DATA.out.amended_json.view()
}

