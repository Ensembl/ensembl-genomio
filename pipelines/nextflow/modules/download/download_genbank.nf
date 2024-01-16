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

process DOWNLOAD_GENBANK {
    tag "${meta.production_name}"
    label 'normal'
    label 'cached'

    input:
        val(meta)

    output:
        tuple val(meta), path("output.gb")

    shell:
    output_file = "output.gb"
    '''
    genbank_download --accession !{meta.accession} --output_file !{output_file} --debug
    '''
    
    stub:
        output_file = "output.gb"
        """
        genbank_download --help
        cp $workflow.projectDir/../../../../data/test/modules/download_genbank/output/*.gb $output_file
        """
}
