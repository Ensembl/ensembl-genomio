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


process COMPARE_GENOME_STATS {
    tag "${db.species}"
    label 'local'

    input:
        tuple val(db), path(ncbi_stats, stageAs: "ncbi_stats.json"),
            path(core_stats, stageAs: "core_stats.json")

    output:
        tuple val(db), val("stats"), path("*_stats.json", includeInputs: true)

    script:
        """
        genome_stats_compare --ncbi $ncbi_stats --core $core_stats --verbose > diff_stats.json
        """
}
