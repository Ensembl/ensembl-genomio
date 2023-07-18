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

process MANIFEST_STATS {
    tag "manifest_stats"
    label 'default'
    publishDir "$out_dir/$accession", mode: 'copy'

    input:
        tuple val(accession), path(genome_files)
        val datasets
        val ignore_failed
        val out_dir

    output:
        tuple val(accession), path("stats.txt")

    script:
        """
        MORE_CMD=""
        if [ "$ignore_failed" = "1" ]; then
            MORE_CMD=" --ignore_failed 1"
        fi
        manifest_stats --manifest_dir "." \
        --datasets_bin "$datasets" \
        --accession "$accession" \
        --stats_file ./stats.txt \
        \$MORE_CMD
        """
}
