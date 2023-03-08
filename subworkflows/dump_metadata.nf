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
include { CHECK_JSON_SCHEMA } from '../modules/check_json_schema.nf'
include { COLLECT_META_FILE } from '../modules/collect_files.nf'
include { PUBLISH_DIR } from '../modules/collect_files.nf'

workflow DUMP_METADATA {
    take:
        server
        db
        filter_map
        out_dir

    emit:
        db

    main:
        seq_regions = DUMP_SEQ_REGIONS(server, db, filter_map, out_dir)
        seq_regions_checked = CHECK_JSON_SCHEMA(seq_regions)

        json_files = seq_regions_checked.merge()
        collect_dir = COLLECT_META_FILE(json_files, db)
        PUBLISH_DIR(collect_dir, db, out_dir)
}
