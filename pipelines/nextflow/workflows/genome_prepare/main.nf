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

if (params.regions_to_exclude) {
    params.regions_to_exclude_list = params.regions_to_exclude?.split(',') as List
}


// Import modules/subworkflows
include { PREPARE_GENOME_METADATA } from '../../modules/prepare_genome_metadata.nf'
include { CHECK_JSON_SCHEMA } from '../../modules/check_json_schema.nf'
include { DOWNLOAD_ASM_DATA } from '../../modules/download_asm_data.nf'
include { UNPACK_FILE } from '../../modules/unpack_gff3.nf'
include { PROCESS_SEQ_REGION } from '../../modules/process_seq_region.nf'


// Run main workflow
workflow {
    ch_metadata_json = PREPARE_GENOME_METADATA(ch_genome_json)
    checked_jsons = CHECK_JSON_SCHEMA(ch_metadata_json)
    checked_jsons.view()
    genome_data = DOWNLOAD_ASM_DATA(checked_jsons)
    ch_process_seqregion = PROCESS_SEQ_REGION(genome_data, params.regions_to_exclude_list)
    ch_process_seqregion.view()
}

