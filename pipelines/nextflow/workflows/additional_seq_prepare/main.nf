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

//default params
params.help = false
params.brc_mode = 1

// mandatory params
params.accession = null
params.prefix = null
params.PROD_NAME = null



// Print usage
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run add_seq_prepare.nf --accession  --prefix "PREFIX_" --PROD_NAME "productionname"

        Mandatory arguments:
        --accession                    A GenBank accession of the sequence you are adding
        --prefix                       Required a string to add to the gene ids, to ensure that they are unique (include  PREFIX_)
        --production_name              Production name of the species

       Optional arguments:
        --output_dir                   Output directory to place final output
        --help                         This usage statement.
        --brc_mode                     By default it is set to 1, set it to 0 if you are not using it for brc
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

assert params.accession, "Parameter 'accession' is not specified"
assert params.prefix, "Parameter 'prefix' is not specified"
assert params.production_name, "Parameter 'production_name' is not specified"

params.meta = [
    "accession": params.accession,
    "production_name": params.production_name,
    "prefix": params.prefix
]

// Import modules/subworkflows
include { additional_seq_prepare } from '../../subworkflows/additional_seq_prepare/main.nf'

// Run main workflow
workflow {
    additional_seq_prepare(params.meta, params.brc_mode, params.output_dir)
}
