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

process DB_FACTORY {
    tag "DB_factory"
    label 'default'

    input:
        val server_map

    output:
        path "*.json"

    script:
        """
        brc_mode=''
        if [ $server_map.brc_mode == 1 ]; then
            brc_mode='--brc_mode 1'
        fi
        db_factory --host '${server_map.host}' --port '${server_map.port}' --user '${server_map.user}' --password '${server_map.password}' --prefix '${server_map.prefix}' \$brc_mode --output_json dbs.json
        """
}
