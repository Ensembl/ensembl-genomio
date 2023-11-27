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

// default params
params.help = false
params.dbname_re = ''
params.output_dir = './gpad_output'
params.password = ''
params.gpad_file = ''

// Print usage
def helpMessage() {
  log.info """
        Mandatory arguments:
        --host, --port, --user         Connection parameters to the SQL servers we getting core db(s) from
        --gpad_file                    Compressed GPAD file to extract data from

        Optional arguments:
        --password                     Password part of the connection parameters
        --dbname_re                    Regexp to match core db name(s) against
        --output_dir                   Name of Output directory to gather prepared outfiles. (default: ${params.output_dir})

        --help                         This usage statement.
        """
}

// Check mandatory parameters
if (params.help) {
    helpMessage()
    exit 0
}
if (!params.gpad_file) {
    exit 0, "Missing gpad_file"
}

def create_server(params) {
    server = [
        "host": params.host,
        "port": params.port,
        "user": params.user,
        "password": ""
    ]
    if (params.password) {
        server["password"] = params.password
    }
    return server
}

def create_filter_map(params) {
    filter_map = [
        "prefix": "",
        "dbname_re": ""
    ]
    if (params.dbname_re) {
        filter_map["dbname_re"] = params.dbname_re
    }
    return filter_map
}

if (params.host && params.port && params.user && params.output_dir) {
    server = create_server(params)
    filter_map = create_filter_map(params)
} else {
    exit 1, "Missing server parameters"
}


include { GPAD_PREPARE } from '../../subworkflows/gpad_prepare_sub/main.nf'

// Run main workflow
workflow {
    GPAD_PREPARE(params.gpad_file, server, filter_map)
}
