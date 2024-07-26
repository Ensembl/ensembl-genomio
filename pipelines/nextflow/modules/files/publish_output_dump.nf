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


process PUBLISH_DIR {
    tag "${db.species}"
    label 'default'
    publishDir "$out_dir/$release_dir/metadata/$db.division/$db.species", mode: 'copy'
    time '5min'

    input:
        tuple val(db), path(data_dir)
        val out_dir
    
    output:
        tuple val(db), path(data_dir, includeInputs: true)

    script:

        // check if the core DB had an expected release version or not (internal/unformatted db name)
        if ( "${db.release}".isEmpty() ) {
            release_dir = "unreleased"
        }
        else{
            release_dir = "build_${db.release}"
        }

        """
        echo "Just copy over the finished files"
        """
}
