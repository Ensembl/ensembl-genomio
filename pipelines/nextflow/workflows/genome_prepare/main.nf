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
include { DATASETS_METADATA } from '../../modules/genome_metadata/datasets_metadata.nf'
include { ACCESSION_METADATA } from '../../modules/genome_metadata/accession_metadata.nf'
// Utilities
include { read_json } from '../../modules/utils/utils.nf'

// function to create meta tuple from genome.json in the form:
//   tuple("accession": accession, "production_name": production_name, "prefix": prefix)
def meta_from_genome_json(json_path) {
    data = read_json(json_path)

    prod_name = data.assembly.accession
    publish_dir = data.assembly.accession
    if ( params.brc_mode ) {
        prod_name = data.BRC4.organism_abbrev
        publish_dir = "${data.BRC4.component}/${data.BRC4.organism_abbrev}"
    } else if ( data.species && data.species.production_name ) {
        prod_name = data.species.production_name
        publish_dir = prod_name
    }
    has_annotation = false
    if (data.annotation) {
        has_annotation = true
    }

    return [
        accession: data.assembly.accession,
        production_name: prod_name,
        publish_dir: publish_dir,
        prefix: "",
        has_annotation: has_annotation,
    ]
}

// Run main workflow
workflow {
    ch_genome_json = Channel.fromPath("${params.input_dir}/*.json", checkIfExists: true)
    accession_meta = ACCESSION_METADATA(ch_genome_json)
    accession_val = accession_meta.map{ accession, meta_file -> accession }
    dataset_report = DATASETS_METADATA(accession_val)
    PREPARE_GENOME_METADATA(accession_meta.join(dataset_report))

    PREPARE_GENOME_METADATA.out.genomic_dataset
        .map{ accession, json_file -> tuple(meta_from_genome_json(json_file), json_file) }
        .set { genome_metadata }

    GENOME_PREPARE(genome_metadata)
}
