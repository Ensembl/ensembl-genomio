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

process DUMP_AGP {
    tag "${db.species}"
    label "variable_2_8_32"
    label "ensembl_scripts_container"
    maxForks params.max_database_forks

    input:
        val db

    output:
        tuple val(db), val("agp"), path("agp_files/*.agp")

    script:
        output_dir = "agp_files"
        password_arg = db.server.password ? "--pass ${db.server.password}" : ""
        """
        # Prepare, make sure we can output something, even if empty
        mkdir $output_dir
        touch $output_dir/empty.agp

        dump_agp.pl \
            --host $db.server.host \
            --port $db.server.port \
            --user $db.server.user \
            $password_arg \
            --dbname $db.server.database \
            -v \
            --output_dir $output_dir

        """
    
    stub:
        output_dir = "agp_files"
        output = "$output_dir/seq.agp"
        """
        mkdir $output_dir
        echo "TEST" > $output
        """
}
