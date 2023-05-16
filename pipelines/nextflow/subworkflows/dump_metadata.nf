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

include { DUMP_SEQ_REGIONS } from '../modules/dump_seq_regions.nf'
include { DUMP_EVENTS } from '../modules/dump_events.nf'
include { CHECK_JSON_SCHEMA } from '../modules/check_json_schema_db.nf'

include { COLLECT_FILES } from '../modules/collect_files_db.nf'
include { MANIFEST } from '../modules/collect_files_db.nf'
include { PUBLISH_DIR } from '../modules/collect_files_db.nf'
include { CHECK_INTEGRITY } from '../modules/check_integrity_db.nf'

workflow DUMP_METADATA {
    take:
        server
        db
        filter_map
        out_dir

    emit:
        db

    main:
        // Generate the files
        seq_regions = DUMP_SEQ_REGIONS(server, db, filter_map)
        seq_regions_checked = CHECK_JSON_SCHEMA(seq_regions)
        events = DUMP_EVENTS(server, db, filter_map)

        // Group the files by db species (use the db object as key)
        // Only keep the files so they are easy to collect
        db_files = seq_regions_checked
            .concat(events)
            .map{ db, name, file_name -> tuple(db, file_name) }
            .groupTuple()

        // Collect, create manifest, and publish
        collect_dir = COLLECT_FILES(db_files)
        manifested_dir = MANIFEST(collect_dir, db)
        manifest_checked = CHECK_INTEGRITY(manifested_dir, filter_map)
        PUBLISH_DIR(manifest_checked, db, out_dir)
}
