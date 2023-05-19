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

// Import modules/subworkflows
include { download_genbank } from '../modules/download_genbank.nf'
include { extract_from_gb } from '../modules/extract_from_gb.nf'
include { PROCESS_GFF3 } from '../modules/process_gff3.nf'
include { GFF3_VALIDATION } from '../modules/gff3_validation.nf'
include { CHECK_JSON_SCHEMA } from '../modules/check_json_schema.nf'
include { JSON_SCHEMA_FACTORY } from '../modules/json_schema_factory.nf'
include { COLLECT_FILES; MANIFEST; PUBLISH_DIR } from '../modules/collect_files.nf'
include { CHECK_INTEGRITY } from '../modules/integrity.nf'
include { MANIFEST_STATS } from '../modules/manifest_stats.nf'

workflow additional_seq_prepare {
    take:
        prefix
        accession
        production_name
        brc_mode
        output_dir
    main:
        // Get the data
        gb_file = download_genbank(accession)

        // Parse data from GB file into GFF3 and json files
        extract_from_gb(prefix, production_name, gb_file)
        PROCESS_GFF3(extract_from_gb.out.gene_gff, extract_from_gb.out.genome, accession)
        GFF3_VALIDATION(PROCESS_GFF3.out.gene_models)

        // Validate files
        json_files = extract_from_gb.out.genome
            .concat(extract_from_gb.out.seq_regions, PROCESS_GFF3.out.functional_annotation)
        CHECK_JSON_SCHEMA(json_files, accession)
        all_files = CHECK_JSON_SCHEMA.out.verified_json
                        .concat(extract_from_gb.out.dna_fasta, extract_from_gb.out.pep_fasta)
                        .map{it -> [accession, it]}
                        .groupTuple()
        
        // Collect in manifest, checks and stats
        collect_dir = COLLECT_FILES(all_files)
        manifest_dired = MANIFEST(collect_dir, accession)
        manifest_checked = CHECK_INTEGRITY(manifest_dired, brc_mode)
        manifest_stated = MANIFEST_STATS(manifest_checked, 'datasets')

        // Publish the data
        PUBLISH_DIR(manifest_stated, output_dir, accession)
}
