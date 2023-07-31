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

process ALLOCATE_IDS {
    label 'local'

    input:
        path(new_registry)
        val(species)
        val(osid)
        path(new_feats)
        val(mode)
    
    output:
        path("output_map.txt")

    script:
    def output_map = "output_map.txt"
    def osid_params = ""
    if (osid.mock) {
        osid_params = "--mock_osid"
    } else {
        exit 1, "Actual OSID not implemented yet: $osid";
    }
    def list = ""
    if (mode == "gene") {
        list = "--gene_list $new_feats"
    } else if (mode == "transcript") {
        list = "--transcript_list $new_feats"
    } else {
        exit 1, "Unsupported allocation mode";
    }
    """
    perl $params.scripts_dir/allocate_stable_ids.pl \\
        --reg ./$new_registry \\
        --species $species \\
        --gene_list $new_feats \\
        --output_map $output_map \\
        $osid_params \\
        --update
    """
}