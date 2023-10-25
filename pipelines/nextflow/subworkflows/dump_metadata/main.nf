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
include { DUMP_EVENTS } from '../../modules/events/dump_events.nf'
include { DUMP_GENOME_META } from '../../modules/genome_metadata/dump_genome_meta.nf'
include { DUMP_GENOME_STATS } from '../../modules/genome_metadata/dump_genome_stats.nf'
include { COMPARE_GENOME_STATS } from '../../modules/genome_metadata/compare_genome_stats.nf'
include { DUMP_NCBI_STATS } from '../../modules/genome_metadata/dump_ncbi_stats.nf'
include { CHECK_JSON_SCHEMA } from '../../modules/schema/check_json_schema_db.nf'

include { COLLECT_FILES } from '../../modules/files/collect_files_db.nf'
include { MANIFEST } from '../../modules/files/collect_files_db.nf'
include { PUBLISH_DIR } from '../../modules/files/collect_files_db.nf'
include { CHECK_INTEGRITY } from '../../modules/manifest/check_integrity_db.nf'

workflow DUMP_METADATA {
    take:
        server
        db
        filter_map

    emit:
        db

    main:
        db_files = Channel.of()

        // Seq regions
        if (params.selection.contains("seq_regions")) {
            seq_regions = DUMP_SEQ_REGIONS(server, db, filter_map)
            seq_regions_checked = CHECK_JSON_SCHEMA(seq_regions)
            db_files = db_files.concat(seq_regions_checked)
        }
        
        // Events
        if (params.selection.contains("events")) {
            events = DUMP_EVENTS(server, db, filter_map)
            db_files = db_files.concat(events)
        }

        // Genome metadata
        if (params.selection.contains("genome_metadata")) {
            genome_meta = DUMP_GENOME_META(server, db, filter_map)
            db_files = db_files.concat(genome_meta)
        }

        // Genome stats
        if (params.selection.contains("stats")) {
            genome_stats = DUMP_GENOME_STATS(server, db)
            ncbi_stats = DUMP_NCBI_STATS(server, db)
            stats = ncbi_stats.join(genome_stats)
            diff_stats = COMPARE_GENOME_STATS(stats)
            stats_files = genome_stats.concat(ncbi_stats).concat(diff_stats)
                .groupTuple()
                .transpose(by: 1)
                .map { db, files -> tuple(db, "stats", files) }
            db_files = db_files.concat(stats_files)
        }

        // Group the files by db species (use the db object as key)
        // Only keep the files so they are easy to collect
        db_files = db_files
            .map{ db, name, file_name -> tuple(db, file_name) }
            .groupTuple()

        // Collect, create manifest, and publish
        collect_dir = COLLECT_FILES(db_files)
        manifested_dir = MANIFEST(collect_dir)
        manifest_checked = CHECK_INTEGRITY(manifested_dir, filter_map)
        PUBLISH_DIR(manifest_checked, params.output_dir)
}
