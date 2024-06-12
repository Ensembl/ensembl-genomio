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


process DOWNLOAD_ASM_DATA {
    tag "${meta.accession}"
    label 'adaptive'
    label 'cached'

    input:
        tuple val(meta), path(json_file)

    output:
        tuple val(meta),
            path("*_assembly_report.txt"),
            path("*_genomic.fna.gz"),
            path("*_genomic.gbff.gz"),
            emit: min_set
        tuple val(meta),
            path("*_genomic.gff.gz"),
            path("*_protein.faa.gz"),
            path("*_genomic.gbff.gz"),
            emit: opt_set, optional: true

    shell:
        '''
        assembly_download --accession !{meta.accession} --download_dir ./ --verbose
        '''
    
    stub:
        """
        assembly_download --help
        cp $workflow.projectDir/../../../../data/test/modules/download_asm_data/output/* .
        """
}
