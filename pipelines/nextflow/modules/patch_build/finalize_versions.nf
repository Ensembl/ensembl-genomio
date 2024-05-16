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

process FINALIZE_VERSIONS {
    label 'local'

    input:
        val(server)
        path(waited_file)
    
    output:
        val("done")

    script:
    """
    function run_sql {
        mysql --host=$server.host --port=$server.port --user=$server.user --password=$server.password \\
        -e "\$1" \\
        $server.database
    }

    # Transfer version from gene down to transcript and translation
    run_sql "UPDATE transcript LEFT JOIN gene USING(gene_id) SET transcript.version = gene.version";
    run_sql "UPDATE translation LEFT JOIN transcript USING(transcript_id) SET translation.version = transcript.version";

    # Set exons to standard name: "{transcript_name}-E{exon_rank}"
    run_sql "UPDATE exon LEFT JOIN exon_transcript USING(exon_id) LEFT JOIN transcript USING(transcript_id) SET exon.stable_id=CONCAT(transcript.stable_id, '-E', rank)"

    # Set exons to standard version 1
    run_sql "UPDATE exon SET version = 1";
    """
}
