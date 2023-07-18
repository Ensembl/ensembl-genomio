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
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_GENOME } from '../../modules/schema/check_json_schema.nf'
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_SEQREG } from '../../modules/schema/check_json_schema.nf'
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_FUNCT } from '../../modules/schema/check_json_schema.nf'
include { DOWNLOAD_ASM_DATA } from '../../modules/download/download_asm_data.nf'
include { UNPACK_GFF3 } from '../../modules/gff3/unpack_gff3.nf'
include { PROCESS_GFF3 } from '../../modules/gff3/process_gff3.nf'
include { GFF3_VALIDATION } from '../../modules/gff3/gff3_validation.nf'
include { PROCESS_SEQ_REGION } from '../../modules/seq_region/process_seq_region.nf'
include { PROCESS_FASTA as PROCESS_FASTA_DNA } from '../../modules/fasta/process_fasta_data.nf'
include { PROCESS_FASTA as PROCESS_FASTA_PEP } from '../../modules/fasta/process_fasta_data.nf'
include { AMEND_GENOME_DATA } from '../../modules/genome_metadata/amend_genome_data.nf'
include { COLLECT_FILES } from '../../modules/files/collect_files.nf'
include { MANIFEST } from '../../modules/manifest/manifest_maker.nf'
include { PUBLISH_DIR } from '../../modules/files/publish_output.nf'
include { CHECK_INTEGRITY } from '../../modules/manifest/integrity.nf'
include { MANIFEST_STATS } from '../../modules/manifest/manifest_stats.nf'


// Run main workflow
workflow GENOME_PREPARE {

    take:
        genomic_dataset // tuple composed of GCA_XXXXXXX.X (as path) and genome.json
        output_dir // User specified or default
        cache_dir
        ignore_failed_stats

    // Main data input to this subworkflow is genomic_dataset tuple
    main:        
        // Verify genome.json schema
        checked_genome = CHECK_JSON_SCHEMA_GENOME(genomic_dataset)

        // Download genome data files. Files may or may not include gene models (gff3) and/or peptides.
        (download_min, download_opt) = DOWNLOAD_ASM_DATA(checked_genome, cache_dir)

        // Check for presence of annotation related output files, then process gene models data if present
        if (download_opt) {

            // Decompress gff3 file, output with accession in tuple
            unpacked_gff = UNPACK_GFF3(download_opt, 'gff')
            
            // Need both gff and genome files in one tuple
            gff_genome = unpacked_gff.concat(checked_genome)
                .groupTuple(size: 2)
                .map{ key, files -> tuple(key, files[0], files[1]) }
            (new_functional_annotation, new_gene_models) = PROCESS_GFF3(gff_genome)

            // Verify functional_annotation.json schema
            functional_annotation = CHECK_JSON_SCHEMA_FUNCT(new_functional_annotation)

            // Tidy and validate gff3 using gff3validator
            gene_models = GFF3_VALIDATION(new_gene_models)

            // Process peptides
            fasta_pep = PROCESS_FASTA_PEP(download_opt, '1')
        }

        // Generate seq_region.json
        new_seq_region = PROCESS_SEQ_REGION(CHECK_JSON_SCHEMA_GENOME.out.verified_json, download_min)

        // Verify seq_region.json schema
        seq_region = CHECK_JSON_SCHEMA_SEQREG(new_seq_region)

        // Process genomic fna
        fasta_dna = PROCESS_FASTA_DNA(download_min, '0')

        // Amend genome data find any additional sequence regions
        amended_genome = AMEND_GENOME_DATA(checked_genome, download_min, params.brc_mode)

        // // Group files
        prepared_files_grouped = amended_genome.concat(
            gene_models, 
            functional_annotation,
            fasta_pep,
            seq_region,
            fasta_dna
        )
        prepared_files = prepared_files_grouped
                    .groupTuple()

        // Collect in manifest, checks and generate sequence stats
        manifest = MANIFEST(prepared_files)
        manifest_checked = CHECK_INTEGRITY(manifest, params.brc_mode)

        // Publish the data to output directory
        PUBLISH_DIR(manifest_checked, output_dir)

        // Include a stats file, outside of the manifest files
        MANIFEST_STATS(manifest_checked, 'datasets', ignore_failed_stats, output_dir)
}
