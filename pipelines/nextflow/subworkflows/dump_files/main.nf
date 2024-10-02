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

include { CHECK_INTEGRITY } from '../../modules/manifest/integrity.nf'
include { COMPARE_GENOME_STATS } from '../../modules/genome_stats/compare_genome_stats.nf'
include { DUMP_AGP } from '../../modules/seq_region/dump_agp.nf'
include { DUMP_ANNOTATION } from '../../modules/annotation/dump_annotation.nf'
include { DUMP_EVENTS } from '../../modules/events/dump_events.nf'
include { DUMP_FASTA_DNA } from '../../modules/fasta/dump_fasta_dna.nf'
include { DUMP_EVENTS } from '../../modules/events/dump_events.nf'
include { DUMP_FASTA_PEPTIDES } from '../../modules/fasta/dump_fasta_peptides.nf'
include { DUMP_GENOME_META } from '../../modules/genome_metadata/dump_genome_meta.nf'
include { DUMP_GENOME_STATS } from '../../modules/genome_stats/dump_genome_stats.nf'
include { DUMP_GFF3 } from '../../modules/gff3/dump_gff3.nf'
include { DOWNLOAD_NCBI_STATS } from '../../modules/download/datasets_genome_meta.nf'
include { DUMP_SEQ_ATTRIB } from '../../modules/seq_region/dump_seq_attrib.nf'
include { DUMP_SEQ_REGIONS } from '../../modules/seq_region/dump_seq_regions.nf'
include { MANIFEST } from '../../modules/manifest/manifest_maker.nf'
include { PUBLISH_DIR } from '../../modules/files/publish_output_dump.nf'

workflow DUMP_FILES {
    take:
        db
        selection
        selection_number

    emit:
        db

    main:
        db_files = Channel.of()

        // Seq regions
        if ("seq_regions" in selection) {
            seq_regions = DUMP_SEQ_REGIONS(db)
            db_files = db_files.mix(seq_regions)
        }

        // AGP
        if ("agp" in selection) {
            agps = DUMP_AGP(db)
            db_files = db_files.mix(agps)
        }

        // Seq region attribs
        if ("seq_attrib" in selection) {
            seq_attrib = DUMP_SEQ_ATTRIB(db)
            db_files = db_files.mix(seq_attrib)
        }

        // Dump DNA sequences
        if ("fasta_dna" in selection) {
            fasta_dna = DUMP_FASTA_DNA(db)
            db_files = db_files.mix(fasta_dna)
        }

        // Dump protein sequences
        if ("fasta_pep" in selection) {
            fasta_pep = DUMP_FASTA_PEPTIDES(db)
            db_files = db_files.mix(fasta_pep)
        }

        // Dump gene models
        if ("gff3" in selection) {
            gff3 = DUMP_GFF3(db)
            db_files = db_files.mix(gff3)
        }

        // Dump functional annotations
        if ("annotation" in selection) {
            annotation = DUMP_ANNOTATION(db)
            db_files = db_files.mix(annotation)
        }
        
        // Events
        if ("events" in selection) {
            events = DUMP_EVENTS(db)
            db_files = db_files.mix(events)
        }

        // Genome metadata
        if ("genome_metadata" in selection) {
            genome_meta = DUMP_GENOME_META(db)
            db_files = db_files.mix(genome_meta)
        }

        // Genome stats
        if ("stats" in selection) {
            genome_stats = DUMP_GENOME_STATS(db)
            ncbi_stats = DOWNLOAD_NCBI_STATS(db)
            stats = ncbi_stats.join(genome_stats)
            stats_files = COMPARE_GENOME_STATS(stats).transpose()
            db_files = db_files.mix(stats_files)
        }

        // Group the files by db species (use the db object as key)
        // Only keep the files so they are easy to collect
        db_files = db_files
            .map{ db, name, file_name -> tuple(db, file_name) }
            .groupTuple(size: selection_number)
            // Flatten in case we get lists of files (e.g. from AGP)
            .map{ db, files -> tuple(db, files.flatten()) }

        // Collect, create manifest, and publish
        manifested_dir = MANIFEST(db_files)
        manifest_checked = CHECK_INTEGRITY(manifested_dir)
        PUBLISH_DIR(manifest_checked, params.output_dir)
}
