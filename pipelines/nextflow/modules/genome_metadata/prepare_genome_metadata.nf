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

// Utilities
include { read_json } from '../../modules/utils/utils.nf'


workflow PREPARE_GENOME_METADATA {
    // Generate a meta value from a db metadata file and an ncbi datasets file
    take:
        input_json  // [ val(meta), path(input_json), path(ncbi_json) ]

    emit:
        genome_metadata
    
    main:
        _PREPARE_GENOME_METADATA(input_json)
            .map{ meta, json_file -> tuple(meta + meta_from_genome_json(json_file), json_file) }
            .set { genome_metadata }
}


process _PREPARE_GENOME_METADATA {
    tag "$meta.id"
    label 'local'

    input:
        tuple val(meta), path(input_json), path(ncbi_json)

    output:
        tuple val(meta), path("genome.json"), emit: genomic_dataset
        
    shell:
    output_json = "genome.json"
    '''
    genome_metadata_prepare --input_file !{input_json} --output_file !{output_json} --ncbi_meta !{ncbi_json}
    '''
}

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
