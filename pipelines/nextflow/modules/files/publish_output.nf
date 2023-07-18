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
    publishDir "$out_dir/$accession", mode: 'copy', overwrite: false
    tag "Publish_${accession}"
    label 'default'
    time '5min'

    input:
        tuple val(accession), path(genome_files)
        val (out_dir)
    
    output:
        path("*.*", includeInputs: true)
    
    script:
        """
        echo "Just copy over the finished files"
        echo "All files here to '$out_dir/$accession'"
        """
}
