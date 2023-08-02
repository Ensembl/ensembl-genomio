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

process CHECK_OSID {
    label 'local'

    input:
        val(osid)
    
    output:
        val(osid)

    script:
    def osid_params = ""
    if (!osid.mock) {
        osid_params = "--osid_url $osid.url --osid_user $osid.user --osid_pass $osid.pass --species $osid.species"
    }
    """
    if [ "$osid_params" != "" ]; then
        perl $params.scripts_dir/osid_check_organisms.pl $osid_params > osid_organism.json
        size=\$(wc -l osid_organism.json)
        if [ \$size -eq 0 ]; then
            echo "Cannot get organism info from OSID for $osid.species"
            exit 1
        fi
    fi
    """
}
