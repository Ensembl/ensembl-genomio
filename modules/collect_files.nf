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

process COLLECT_META_FILE {
    publishDir "$out_dir/metadata/$db.division/$db.species", mode: 'copy'
    tag "Collect_file_${name}"
    label 'default'
    time '5min'

    input:
        tuple val(name), path(meta_file)
        val db
        val out_dir
    
    output:
        path meta_file
    
    script:
        """
        echo "Copy final file ${meta_file}"
        """
}
