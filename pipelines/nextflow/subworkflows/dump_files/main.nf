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

include { DUMP_SEQ_REGIONS } from '../../modules/seq_region/dump_seq_regions.nf'
include { DUMP_FASTA_DNA } from '../../modules/fasta/dump_fasta_dna.nf'
include { DUMP_FASTA_PEPTIDES } from '../../modules/fasta/dump_fasta_peptides.nf'
// include { DUMP_GFF3 } from '../../modules/gff3/dump_gff3.nf'
include { DUMP_EVENTS } from '../../modules/events/dump_events.nf'
include { DUMP_GENOME_META } from '../../modules/genome_metadata/dump_genome_meta.nf'
include { DUMP_GENOME_STATS } from '../../modules/genome_stats/dump_genome_stats.nf'
include { COMPARE_GENOME_STATS } from '../../modules/genome_stats/compare_genome_stats.nf'
include { DUMP_NCBI_STATS } from '../../modules/genome_metadata/dump_ncbi_stats.nf'

include { MANIFEST } from '../../modules/manifest/manifest_maker.nf'
include { CHECK_INTEGRITY } from '../../modules/manifest/integrity.nf'
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
            db_files = db_files.concat(seq_regions)
        }

        // Dump DNA sequences
        if ("fasta_dna" in selection) {
            fasta_dna = DUMP_FASTA_DNA(db)
            db_files = db_files.concat(fasta_dna)
        }

        // Dump protein sequences
        if ("fasta_pep" in selection) {
            fasta_pep = DUMP_FASTA_PEPTIDES(db)
            db_files = db_files.concat(fasta_pep)
        }

        // Dump gene models
        // if ("gff3" in selection) {
            // gff3 = DUMP_GFF3(db)
            // db_files = db_files.concat(gff3)
        // }
        
        // Events
        if ("events" in selection) {
            events = DUMP_EVENTS(db)
            db_files = db_files.concat(events)
        }

        // Genome metadata
        if ("genome_metadata" in selection) {
            genome_meta = DUMP_GENOME_META(db)
            db_files = db_files.concat(genome_meta)
        }

        // Genome stats
        if ("stats" in selection) {
            genome_stats = DUMP_GENOME_STATS(db)
            ncbi_stats = DUMP_NCBI_STATS(db)
            stats = ncbi_stats.join(genome_stats)
            stats_files = COMPARE_GENOME_STATS(stats).transpose()
            db_files = db_files.concat(stats_files)
        }

        // Group the files by db species (use the db object as key)
        // Only keep the files so they are easy to collect
        db_files = db_files
            .map{ db, name, file_name -> tuple(db, file_name) }
            .groupTuple(size: selection_number)

        // Collect, create manifest, and publish
        manifested_dir = MANIFEST(db_files)
        manifest_checked = CHECK_INTEGRITY(manifested_dir)
        PUBLISH_DIR(manifest_checked, params.output_dir)
}
