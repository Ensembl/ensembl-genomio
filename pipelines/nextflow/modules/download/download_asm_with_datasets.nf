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


process DOWNLOAD_ASM_NCBI_DATSETS {
    tag "${meta.accession}"
    label 'adaptive'
    label 'cached'
    label 'datasets_cli'

    input:
        tuple val(meta), path(json_file)

    output:
        tuple val(meta),
            path("*_assembly_report.txt"),
            path("*_genomic.fna.gz"),
            path("genomic.gbff.gz"),
            emit: min_set
        tuple val(meta),
            path("genomic.gff.gz"),
            path("protein.faa.gz"),
            path("genomic.gbff.gz"),
            emit: opt_set, optional: true

    shell:
        '''
        datasets download genome accession !{meta.accession} \
        --include seq-report,genome,gbff,gff3,protein \
        --assembly-source all \
        --dehydrated \
        --filename !{meta.accession}_ncbi.zip \
        --no-progressbar

        unzip !{meta.accession}_ncbi.zip -d !{meta.accession}_ncbi

        datasets rehydrate --directory !{meta.accession}_ncbi \
        --gzip \
        --no-progressbar

        cp -f !{meta.accession}_ncbi/ncbi_dataset/data/!{meta.accession}/* ./

        dataformat tsv genome-seq --inputfile sequence_report.jsonl > !{meta.accession}_assembly_report.txt
        '''

    stub:
        """
        datasets --help
        cp $workflow.projectDir/../../../../data/test/modules/download_asm_data/output/* .
        """
}