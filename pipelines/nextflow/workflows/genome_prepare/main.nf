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

// Check mandatory parameters
if (params.input_dir) {
    ch_genome_json = Channel.fromPath("${params.input_dir}/*.json", checkIfExists: true)
} else {
    exit 1, 'Input directory not specified!'
}


// Import modules/subworkflows
include { PREPARE_GENOME_METADATA } from '../../modules/prepare_genome_metadata.nf'


// Run main workflow
workflow {
    ch_metadata_json = PREPARE_GENOME_METADATA(ch_genome_json)
}