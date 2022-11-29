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

// Import modules/subworkflows
include { JSON_SCHEMA_FACTORY } from '../modules/json_schema_factory.nf'
include { CHECK_JSON_SCHEMA } from '../modules/check_json_schema.nf'


// Import utility functions
include { get_key_list } from '../modules/utils.nf'


workflow MY_SUBWORKFLOW {
    take:
        ch_manifest_dir

    // emit:
    //     Nothing to emit

    main:
        metadata_types = get_key_list(params.json_schemas)
        JSON_SCHEMA_FACTORY(ch_manifest_dir, metadata_types)
        CHECK_JSON_SCHEMA(JSON_SCHEMA_FACTORY.out.flatten())
}
