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
        val(finalized)
    
    output:
        path("patch_ids_error_report.log")

    shell:
    '''
    error_report="patch_ids_error_report.log"
    touch $error_report

    function run_sql {
        mysql --host=!{server.host} \\
        --port=!{server.port} \\
        --user=!{server.user} \\
        --password=!{server.password} \\
        -e "$1" \\
        !{server.database}
    }
    
    function check_duplicates {
        table=$1
        dups=$(run_sql "SELECT stable_id FROM $table GROUP BY stable_id HAVING count(*)>1" | wc -l)
        if [ "$dups" -gt "0" ]; then
            echo "$dups duplicate $table stable ids" >> $error_report
            echo 1
            return
        fi
        echo 0
    }
    
    function check_apollo_id {
        table=$1
        apollo_pattern="[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
        apids=$(run_sql "SELECT stable_id FROM $table" | grep -E "$apollo_pattern" | wc -l)
        if [ "$apids" -gt "0" ]; then
            echo "$apids Apollo IDs remain in $table" >> $error_report
            echo 1
            return
        fi
        echo 0
    }

    # Check for stable ids duplicates
    errors=0
    for table in "gene" "transcript" "translation"; do
        error=$(check_duplicates $table)
        errors=$(($errors + $error))
    done

    # Check for Apollo ids remaining
    for table in "gene" "transcript" "translation" "exon"; do
        error=$(check_apollo_id $table)
        errors=$(($errors + $error))
    done

    if [ "$errors" -gt "0" ]; then
        echo "$errors errors found, abort" >> $error_report
        cat $error_report
        exit 1
    fi
    '''
}
