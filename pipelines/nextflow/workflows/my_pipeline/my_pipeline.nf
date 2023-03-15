#!/usr/bin/env nextflow
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

// Check mandatory parameters
if (params.manifest_dir) {
    ch_manifest_dir = file(params.manifest_dir, checkIfExists: true)  // equivalent to: Channel.fromPath(...)
} else {
    exit 1, 'Manifest directory not specified!'
}


// Import modules/subworkflows
include { MY_SUBWORKFLOW } from '../../subworkflows/my_subworkflow.nf'


// Run main workflow
workflow {
    MY_SUBWORKFLOW(ch_manifest_dir)
}
