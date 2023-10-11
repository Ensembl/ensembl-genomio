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
include { DOWNLOAD_GENBANK } from '../../modules/download/download_genbank.nf'
include { EXTRACT_FROM_GB } from '../../modules/genbank/extract_from_gb.nf'
include { PROCESS_GFF3 } from '../../modules/gff3/process_gff3.nf'
include { GFF3_VALIDATION } from '../../modules/gff3/gff3_validation.nf'
include { CHECK_JSON_SCHEMA } from '../../modules/schema/check_json_schema.nf'
include { MANIFEST } from '../../modules/manifest/manifest_maker.nf'
include { PUBLISH_DIR } from '../../modules/files/publish_output.nf'
include { CHECK_INTEGRITY } from '../../modules/manifest/integrity.nf'
include { MANIFEST_STATS } from '../../modules/manifest/manifest_stats.nf'

workflow additional_seq_prepare {
    take:
        meta
        brc_mode
        output_dir
    main:
        // We expect every input and output stream to have `meta` as the first val in the form of:
        //   tuple("accession": accession, "production_name": production_name, "prefix": prefix)

        // Get the data
        gb_file = DOWNLOAD_GENBANK(meta)

        // Parse data from GB file into GFF3 and json files
        (gff_genome, gb_genome, gb_seq_regions, gb_dna_fasta, gb_pep_fasta) = EXTRACT_FROM_GB(gb_file)

        // Process the GFF and GB files into a cleaned GFF and a functional_annotation file
        (new_functional_annotation, new_gene_models) = PROCESS_GFF3(gff_genome.join(gb_genome))

        // // Tidy and validate gff3 using gff3validator
        gene_models = GFF3_VALIDATION(new_gene_models)

        // Validate files
        json_files = gb_genome.concat(
            gb_seq_regions,
            new_functional_annotation
        )
        
        // Verify schemas 
        checked_json_files = CHECK_JSON_SCHEMA(json_files)

        // Gather json and fasta files and reduce to unique input accession
        all_files = checked_json_files.concat(
            gene_models,
            gb_dna_fasta,
            gb_pep_fasta
        ).groupTuple()
        
        manifest_bundle = MANIFEST(all_files)
        
        // Checks and generate sequence stats for manifest
        manifest_checked = CHECK_INTEGRITY(manifest_bundle, params.brc_mode)
        
        manifest_stated = MANIFEST_STATS(manifest_checked, 'datasets', 0)

        // Publish the data to output directory
        PUBLISH_DIR(manifest_stated, output_dir)
}
