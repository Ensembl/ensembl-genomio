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

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

if (params.help) {
     log.info paramsHelp("nextflow run additional_seq_prepare/main.nf --accession 'GB_ACCESSION' --prefix 'PREFIX_' --production_name 'SPECIES_PROD_NAME'")
    exit 0
}

validateParameters()
log.info paramsSummaryLog(workflow)

if (params.brc_mode) {
    params.brc_mode = params.brc_mode as Integer
}

meta = [
    "accession": params.accession,
    "production_name": params.production_name,
    "prefix": params.prefix,
    "publish_dir": params.accession,
]

// Import modules/subworkflows
include { additional_seq_prepare } from '../../subworkflows/additional_seq_prepare/main.nf'

// Run main workflow
workflow {
    additional_seq_prepare(meta)
}
