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

// Import modules/subworkflows
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_GENOME } from '../../modules/check_json_schema.nf'
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_SEQREG } from '../../modules/check_json_schema.nf'
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_FUNCT } from '../../modules/check_json_schema.nf'
include { DOWNLOAD_ASM_DATA } from '../../modules/download_asm_data.nf'
include { UNPACK_GFF3 } from '../../modules/unpack_gff3.nf'
include { PROCESS_GFF3 } from '../../modules/process_gff3.nf'
include { GFF3_VALIDATION } from '../../modules/gff3_validation.nf'
include { PROCESS_SEQ_REGION } from '../../modules/process_seq_region.nf'
include { PROCESS_FASTA as PROCESS_FASTA_DNA } from '../../modules/process_fasta_data.nf'
include { PROCESS_FASTA as PROCESS_FASTA_PEP } from '../../modules/process_fasta_data.nf'
include { AMEND_GENOME_DATA } from '../../modules/amend_genome_data.nf'
include { COLLECT_FILES } from '../../modules/collect_files.nf'
include { MANIFEST } from '../../modules/manifest_maker.nf'
include { PUBLISH_DIR } from '../../modules/publish_output.nf'
include { CHECK_INTEGRITY } from '../../modules/integrity.nf'
include { MANIFEST_STATS } from '../../modules/manifest_stats.nf'


// Run main workflow
workflow GENOME_PREPARE {

    take:
        genomic_dataset // tuple composed of GCA_XXXXXXX.X (as path) and genome.json
        output_dir // User specified or default
        cache_dir

    // Main data input to this subworkflow is genomic_dataset tuple
    main:        
        // Verify genome.json schema
        CHECK_JSON_SCHEMA_GENOME(genomic_dataset)

        // Download genome data files. Files may or may not include gene models (gff3) and/or peptides.
        DOWNLOAD_ASM_DATA(CHECK_JSON_SCHEMA_GENOME.out.verified_json, cache_dir)

        // Check for presence of annotation related output files, then process gene models data if present
        if (DOWNLOAD_ASM_DATA.out.opt_set) {

            // Decompress gff3 file, output with accession in tuple
            unpacked_gff = UNPACK_GFF3(DOWNLOAD_ASM_DATA.out.opt_set, 'gff')
            
            PROCESS_GFF3(unpacked_gff, CHECK_JSON_SCHEMA_GENOME.out.verified_json)

            // Verify functional_annotation.json schema
            CHECK_JSON_SCHEMA_FUNCT(PROCESS_GFF3.out.functional_annotation)

            // Tidy and validate gff3 using gff3validator (NOTE: Requires `module load libffi-3.3-gcc-9.3.0-cgokng6`)
            GFF3_VALIDATION(PROCESS_GFF3.out.functional_annotation, PROCESS_GFF3.out.gene_models)

            // Process peptides
            PROCESS_FASTA_PEP(DOWNLOAD_ASM_DATA.out.opt_set, '1')
        }

        // Generate seq_region.json
        PROCESS_SEQ_REGION(CHECK_JSON_SCHEMA_GENOME.out.verified_json, DOWNLOAD_ASM_DATA.out.min_set)

        // Verify seq_region.json schema
        CHECK_JSON_SCHEMA_SEQREG(PROCESS_SEQ_REGION.out.seq_region)

        // Process genomic fna
        PROCESS_FASTA_DNA(DOWNLOAD_ASM_DATA.out.min_set, '0')

        // Amend genome data find any additional sequence regions
        AMEND_GENOME_DATA(CHECK_JSON_SCHEMA_GENOME.out.verified_json, DOWNLOAD_ASM_DATA.out.min_set, params.brc_mode)

        // // Group files
        prepared_files = AMEND_GENOME_DATA.out.amended_json
                    .concat(PROCESS_GFF3.out.gene_models, 
                    CHECK_JSON_SCHEMA_FUNCT.out.verified_json,
                    PROCESS_FASTA_PEP.out.processed_fasta, 
                    CHECK_JSON_SCHEMA_SEQREG.out.verified_json,
                    PROCESS_FASTA_DNA.out.processed_fasta)
                    .groupTuple()

        // Collect in manifest, checks and generate sequence stats
        collect_dir = COLLECT_FILES(prepared_files)

        manifest_dired = MANIFEST(collect_dir)
        
        manifest_checked = CHECK_INTEGRITY(manifest_dired, params.brc_mode)
        
        manifest_stated = MANIFEST_STATS(manifest_checked, 'datasets')

        // Publish the data to output directory
        PUBLISH_DIR(manifest_stated, output_dir)
}