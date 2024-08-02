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

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
if (params.help) {
     log.info paramsHelp("nextflow run genome_prepare/main.nf --input_dir <json_input_dir>")
    exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

if (params.brc_mode) {
    params.brc_mode = params.brc_mode as Integer
}

// Import subworkflow
include { GENOME_PREPARE } from '../../subworkflows/genome_prepare/main.nf'
// Import module
include { PREPARE_GENOME_METADATA } from '../../modules/genome_metadata/prepare_genome_metadata.nf'
include { DOWNLOAD_NCBI_STATS } from '../../modules/download/datasets_genome_meta.nf'
include { ACCESSION_METADATA } from '../../modules/genome_metadata/accession_metadata.nf'

// Run main workflow
workflow {
    ch_genome_json = Channel.fromPath("${params.input_dir}/*.json", checkIfExists: true)
    accession_meta = ACCESSION_METADATA(ch_genome_json)
    dataset_report = DOWNLOAD_NCBI_STATS(accession_meta.map{ meta, json -> meta })
    genome_metadata = PREPARE_GENOME_METADATA(accession_meta.join(dataset_report))

    GENOME_PREPARE(genome_metadata)
}
