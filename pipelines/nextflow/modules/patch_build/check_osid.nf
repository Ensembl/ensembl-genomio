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
    def osid_file = "osid_organism.json"
    """
    if [ "$osid.mock" != "true" ]; then
        python $params.scripts_dir/osid_check_organisms.py \\
        --url $osid.url \\
        --user $osid.user \\
        --key $osid.pass \\
        --species $osid.species \\
         > $osid_file
        FOUND=\$(cat osid_organism.json)
        if [ "\$FOUND" = "[]" ]; then
            echo "Cannot get organism info from OSID for $osid.species"
            exit 1
        else
            echo "$osid.species is in OSID"
        fi
    fi
    """
}
