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

process TRANSFER_IDS {
    label 'local'

    input:
        path(changed_genes)
        path(old_registry, stageAs: "old_registry.pm")
        path(new_registry, stageAs: "new_registry.pm")
        val(species)

    output:
        path("new_transcripts.txt"), emit: new_transcripts
        path("mapped_transcripts.txt"), emit: mapped_transcripts

    script:
    def new_transcripts = "new_transcripts.txt"
    def mapped_transcripts = "mapped_transcripts.txt"
    """
    perl $params.scripts_dir/transfer_ids.pl \\
        --mapping $changed_genes \\
        --old ./$old_registry \\
        --new ./$new_registry \\
        --species $species \\
        --mapped_transcripts $mapped_transcripts \\
        --missed_transcripts $new_transcripts \\
        --update
    """
}
