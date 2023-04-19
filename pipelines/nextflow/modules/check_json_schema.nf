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

process CHECK_JSON_SCHEMA {
    label 'default'

    input:
    val metadata_type
    path accession_dir

    output:
    // tuple val(accession), path("${accession_dir}/genome.json")
    val "$accession", emit: gca
    path "${accession_dir}/${metadata_type}.json", emit: json_file

    script:
    accession = accession_dir.getName()
    json_schema = params.json_schemas["$metadata_type"]
    """
    check_json_schema --json_file ${accession_dir}/${metadata_type}.json --json_schema $json_schema
    """
}
