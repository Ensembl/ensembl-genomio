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
    tag "${db.species}"
    label "variable_2_8_32"

    input:
        val server
        val db
        val filter_map

    output:
        tuple val(db), val("seq_region"), path("seq_region.json")

    script:
        """
        brc_mode=''
        if [ $filter_map.brc_mode == 1 ]; then
            brc_mode='--brc_mode 1'
        fi
        touch seq_region.json
        seq_region_dumper --host '${server.host}' \
            --port '${server.port}' \
            --user '${server.user}' \
            --password '${server.password}' \
            --database '${db.database}' \
            \$brc_mode \
            --output_json seq_region.json
        """
}
