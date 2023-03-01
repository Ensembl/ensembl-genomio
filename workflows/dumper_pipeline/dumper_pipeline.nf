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
        "brc_mode": 0,
        "prefix": ""
    ]
    if (params.brc_mode) {
        filter_map["brc_mode"] = 1
    }
    if (params.prefix) {
        filter_map["prefix"] = params.prefix
    }
    return filter_map
}

if (params.host && params.port && params.user && params.out_dir) {
    server = create_server(params)
    filter_map = create_filter_map(params)
} else {
    exit 1, "Missing server parameters"
}

include { DUMP_SQL } from '../../subworkflows/dump_sql.nf'

// Run main workflow
workflow {
    DUMP_SQL(server, filter_map, params.out_dir)
}
