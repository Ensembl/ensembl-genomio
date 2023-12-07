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

process SCAN_GPAD {
    label "long"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        path gpad_file
        path uniprot_maps
    
    output:
        path("*.gpad")

    script:
        all_maps = "all_maps.tsv"
        """
        cat $uniprot_maps > $all_maps
        xrefs_scan_gpad --map $all_maps --gpad $gpad_file --output_dir ./ -v
        """
    
    stub:
        all_maps = "all_maps.tsv"
        output = "mapped.gpad"
        """
        touch species_species.gpad
        cat $uniprot_maps > $all_maps
        touch $output
        """
}
