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

include { DB_FACTORY } from '../modules/db_factory.nf'
include { DUMP_DB } from '../modules/dump_db.nf'
include { read_json } from '../modules/utils.nf'

workflow DUMP_SQL {
    take:
        server
        filter_map
        out_dir

    emit:
        dbs

    main:
        dbs = DB_FACTORY(server, filter_map)
            .map(it -> read_json(it))
            .flatten()
        
        DUMP_DB(server, dbs, out_dir)
}
