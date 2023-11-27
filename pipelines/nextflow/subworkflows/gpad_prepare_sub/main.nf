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

include { DB_FACTORY } from '../../modules/database/db_factory.nf'
include { read_json } from '../../modules/utils/utils.nf'
include { DUMP_UNIPROT_MAP } from '../../modules/xrefs/dump_uniprot_map.nf'
include { SCAN_GPAD } from '../../modules/xrefs/scan_gpad.nf'

workflow GPAD_PREPARE {
    take:
        gpad_file
        server
        filter_map

    main:
        db = DB_FACTORY(server, filter_map)
            .map(it -> read_json(it))
            .flatten()

    uniprot = DUMP_UNIPROT_MAP(db)
    uniprots = uniprot.collect().transpose()
    uniprot_map = SCAN_GPAD(gpad_file, uniprots)
}
