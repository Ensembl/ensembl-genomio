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

process CHECK_PATCH {
    label 'local'

    input:
        val(server)
        path(waited_file)

    shell:
        '''
        function mysql_command {
            mysql !{server.database} \
                --host !{server.host} \
                --port !{server.port} \
                --user !{server.user} \
                --password !{server.password} \
                -e "$1"
        }

        # Purge all current canonical transcripts
        mysql_command "UPDATE gene SET canonical_transcript_id = 0"
        mysql_command "DELETE transcript_attrib FROM transcript_attrib LEFT JOIN attrib_type USING(attrib_type_id) WHERE code='is_canonical'"

        perl $ENSEMBL_ROOT_DIR/ensembl/misc-scripts/canonical_transcripts/select_canonical_transcripts.pl \
            --dbhost !{server.host} \\
            --dbport !{server.port} \\
            --dbuser !{server.user} \\
            --dbpassword !{server.password} \\
            --dbname !{server.database} \\
            -coord_system toplevel \
            -write
        '''
}
