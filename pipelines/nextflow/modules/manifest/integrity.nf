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
    tag "${accession}"
    label 'variable_2_8_32'
    errorStrategy 'finish'
    time '1h'

    input:
        tuple val(accession), path(genome_files)
        val brc_mode
    
    output:
        tuple val(accession), path("*.*", includeInputs: true)
    
    script:
        """
        brc_mode=''
        if [ $brc_mode == 1 ]; then
            brc_mode='--brc_mode 1'
        fi
        check_integrity --manifest_file ./manifest.json \
            \$brc_mode
        """
}
