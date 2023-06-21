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

// dealing with custom $NXF_WORK based workDir
//   don't do anything if "-work-dir (-w)" option specified on command line
cmd_line = binding.variables.workflow.commandLine
cmd_line_has_wd = cmd_line.contains(" -w ") || cmd_line.contains(" -work_dir ")
if (!cmd_line_has_wd) {
  log.info " no -work-dir (-w) option specified. Trying to build one base on NXF_WORK env"
  nxf_work_env = binding.getVariable('NXF_WORK')
  if (nxf_work_env) {
    session.workDir = ("${nxf_work_env}/dumper_pipeline" as Path).complete()
    session.workDir.mkdirs()
    workDir = session.workDir as String
  }
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
        "brc_mode": 0,
        "prefix": ""
    ]
    if (params.brc_mode) {
        filter_map["brc_mode"] = 1
    }
    if (params.prefix) {
        filter_map["prefix"] = params.prefix
    }
    if (params.dbname_re) {
        filter_map["dbname_re"] = params.dbname_re
    }
    return filter_map
}

if (params.host && params.port && params.user && params.out_dir) {
    server = create_server(params)
    filter_map = create_filter_map(params)
} else {
    exit 1, "Missing server parameters"
}

include { DUMP_SQL } from '../../subworkflows/dump_sql/main.nf'
include { DUMP_METADATA } from '../../subworkflows/dump_metadata/main.nf'
include { DB_FACTORY } from '../../modules/database/db_factory.nf'
include { read_json } from '../../modules/utils/utils.nf'

// Run main workflow
workflow {
    dbs = DB_FACTORY(server, filter_map)
        .map(it -> read_json(it))
        .flatten()
    DUMP_SQL(server, dbs, filter_map, params.out_dir)
    DUMP_METADATA(server, dbs, filter_map, params.out_dir)
}
