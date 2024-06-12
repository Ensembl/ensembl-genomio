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

process DUMP_GFF3 {
    tag "${db.species}"
    label "variable_2_8_32"
    maxForks params.max_database_forks

    input:
        val db

    output:
        tuple val(db), val("gff3"), path("*.gff3")

    script:
        output = "${db.species}.gff3"
        registry = "registry.pm"
        """
        # Create a registry
        echo "use strict;
        use warnings;
        use Bio::EnsEMBL::Registry;

        {
        Bio::EnsEMBL::Registry->load_registry_from_db(
            -host       => '$db.server.host',
            -port       => $db.server.port,
            -user       => '$db.server.user',
            -pass       => '$db.server.password',
            -db_name  => '$db.server.database',
            -species  => '$db.production_name',
        );
        }

        1;" > $registry

        standaloneJob.pl Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile -reg_conf './$registry' -input_id "{
            species =>           '$db.production_name',
            feature_type => ['Gene', 'Transcript', 'Translation'],
            per_chromosome =>    0,
            include_scaffold =>  1,
            db_type =>           'core',
            xrefs =>             0,
            release => $db.release,
            base_path => '.'
        }"
        mv *.gff3 $output
        """
    
    stub:
        output = "gene_models.gff3"
        """
        echo "##gff-version 3\n" > $output
        """
}
