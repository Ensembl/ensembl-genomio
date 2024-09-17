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
include { AMEND_GENOME_DATA } from '../../modules/genome_metadata/amend_genome_data.nf'
include { CHECK_INTEGRITY } from '../../modules/manifest/integrity.nf'
include { CHECK_JSON_SCHEMA as CHECK_JSON_SCHEMA_GENOME } from '../../modules/schema/check_json_schema.nf'
include { DOWNLOAD_ASM_DATA } from '../../modules/download/download_asm_data.nf'
include { GFF3_VALIDATION } from '../../modules/gff3/gff3_validation.nf'
include { MANIFEST } from '../../modules/manifest/manifest_maker.nf'
include { MANIFEST_STATS } from '../../modules/manifest/manifest_stats.nf'
include { PROCESS_FASTA as PROCESS_FASTA_DNA } from '../../modules/fasta/process_fasta_data.nf'
include { PROCESS_FASTA as PROCESS_FASTA_PEP } from '../../modules/fasta/process_fasta_data.nf'
include { PROCESS_GFF3 } from '../../modules/gff3/process_gff3.nf'
include { PROCESS_SEQ_REGION } from '../../modules/seq_region/process_seq_region.nf'
include { PUBLISH_DIR } from '../../modules/files/publish_output.nf'
include { UNPACK_GFF3 } from '../../modules/gff3/unpack_gff3.nf'


// Main workflow
workflow GENOME_PREPARE {
    // We expect every input and output stream to have `meta` as the first element in the form of:
    //   tuple("accession": accession, "production_name": production_name, "prefix": prefix)

    take:
        genomic_dataset // tuple(`meta`, path("genome.json"))

    main:
        // Verify genome.json schema
        checked_genome = CHECK_JSON_SCHEMA_GENOME(genomic_dataset).verified_json

        // Download genome data files. Files may or may not include gene models (GFF3) and/or peptides.
        DOWNLOAD_ASM_DATA(checked_genome)
        download_min = DOWNLOAD_ASM_DATA.out.min_set
        download_opt = DOWNLOAD_ASM_DATA.out.opt_set

        // Decompress GFF3 file, output with accession in tuple
        unpacked_gff = UNPACK_GFF3(download_opt, 'gff')

        // Process the GB and GFF3 files into a cleaned GFF3 and a functional_annotation files
        genome_gff_files = checked_genome.join(unpacked_gff, failOnDuplicate: true)
        PROCESS_GFF3(genome_gff_files)
        functional_annotation = PROCESS_GFF3.out.functional_annotation
        new_gene_models = PROCESS_GFF3.out.gene_models

        // Tidy and validate gff3 using gff3validator
        gene_models = GFF3_VALIDATION(new_gene_models).gene_models

        // Process peptides
        fasta_pep = PROCESS_FASTA_PEP(download_opt, 1).processed_fasta

        // Group all the genome data files under the same meta key
        genome_data_files = checked_genome.join(download_min, failOnDuplicate: true, failOnMismatch: true)

        // Generate seq_region.json
        seq_region = PROCESS_SEQ_REGION(genome_data_files).seq_region

        // Process genomic fna
        fasta_dna = PROCESS_FASTA_DNA(download_min, 0).processed_fasta

        // Amend genome data find any additional sequence regions
        amended_genome = AMEND_GENOME_DATA(genome_data_files).amended_json

        // Group files
        prepared_files = amended_genome.mix(
                seq_region,
                fasta_dna,
                gene_models,
                functional_annotation,
                fasta_pep,
            )
            .map{ meta, file ->
                key = meta
                if (meta["has_annotation"]) {
                    key = groupKey(meta, 6)
                } else {
                    key = groupKey(meta, 3)
                }
                [key, file] }
            .groupTuple()

        // Checks and generate sequence stats for manifest
        manifest_bundle = MANIFEST(prepared_files)

        manifest_checked = CHECK_INTEGRITY(manifest_bundle)
        
        manifest_stated = MANIFEST_STATS(manifest_checked, 'datasets', params.ncbi_check)

        // Publish the data to output directory
        PUBLISH_DIR(manifest_stated, params.output_dir)
}
