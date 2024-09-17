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


process CHECK_INTEGRITY {
    tag "${meta.accession}"
    label 'variable_2_8_32'
    time '1h'

    input:
        tuple val(meta), path(manifest_files)
    
    output:
        tuple val(meta), path("*.*", includeInputs: true)

    shell:
        integrity_file = "integrity.out"
        '''
        manifest_check_integrity \
            --manifest_file ./manifest.json \
            --no_fail \
            > !{integrity_file}
        
        # Only keep integrity file if there are errors to report
        if [ ! -s !{integrity_file} ]
            then rm !{integrity_file}
        fi
        '''
}
