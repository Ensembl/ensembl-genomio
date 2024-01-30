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
    tag "publish_${meta.accession}"
    label 'default'
    publishDir "$out_dir/${meta.publish_dir}", mode: 'copy', overwrite: false
    time '5min'

    input:
        tuple val(meta), path(bundle_files)
        val (out_dir)
    
    output:
        path "*", includeInputs: true
    
    shell:
        '''
        echo "Just publishing the finished files"
        echo "To '!{out_dir}/!{meta.accession}' for accession '!{meta.accession}'"
        '''
}
