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

process DUMP_SEQ_REGIONS {
    tag "Dump_seq_regions:${db.species}"
    label 'default'
    time '1h'

    input:
        val server
        val db
        val filter_map
        val out_dir

    output:
        path "seq_region.json"

    script:
        """
        brc_mode=''
        if [ $filter_map.brc_mode == 1 ]; then
            brc_mode='--brc_mode 1'
        fi
        seq_region_dumper --host '${server.host}' \
            --port '${server.port}' \
            --user '${server.user}' \
            --password '${server.password}' \
            --database '${db.database}' \
            \$brc_mode \
            --output_json seq_region.json
        """
}
