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

// default params
params.help = false

// Print usage
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the 'Genome prepare' pipeline is as follows:

        nextflow run workflows/genome_prepare/main.nf --input_dir <json_input_dir>

        Mandatory arguments:
        --input_dir                    Location of input json(s) with component/organism genome metadata
        --cache_dir                    Cache dir for the downloaded assembly data

        Optional arguments:
        --output_dir                   Name of Output directory to gather prepared outfiles. Default -> 'Output_GenomePrepare'.
        --ignore_failed_stats          Skip and continue if stats counts fail
        --help                         This usage statement.
        """
}

// Check mandatory parameters
if (params.help) {
    helpMessage()
    exit 0
}
if (params.input_dir) {
    ch_genome_json = Channel.fromPath("${params.input_dir}/*.json", checkIfExists: true)
} else {
    exit 1, 'Input directory not specified!'
}
if (!params.cache_dir) {
    params.cache_dir = "./genome_prepare_download_cache"
}
if (!params.ignore_failed_stats) {
    params.ignore_failed_stats = 1
}

// Import subworkflow
include { GENOME_PREPARE } from '../../subworkflows/genome_prepare/main.nf'
// Import module
include { PREPARE_GENOME_METADATA } from '../../modules/genome_metadata/prepare_genome_metadata.nf'


// Run main workflow
workflow {
    PREPARE_GENOME_METADATA(ch_genome_json)
    PREPARE_GENOME_METADATA.out.genomic_dataset
        .map{ gca_dir, json_file -> tuple( gca_dir.getName(), json_file ) }
        .set { genome_metadata }
    GENOME_PREPARE(genome_metadata, params.output_dir, params.cache_dir, params.ignore_failed_stats)
}
