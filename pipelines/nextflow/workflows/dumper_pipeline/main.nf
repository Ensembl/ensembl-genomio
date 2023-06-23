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
include { updateWorkDirAsNeeded } from '../../modules/utils/utils.nf'
updateWorkDirAsNeeded("dumper_pipeline2")


// default params
params.help = false
params.prefix = ''
params.brc_mode = 0
params.dbname_re = ''
params.output_dir = './dumper_output'
params.password = ''

// Print usage
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the 'Dumper' pipeline is as follows:

        CMD=<dba_alias>
        pushd data
          data_dir=\$(pwd)
          nextflow run \\
            -w \${data_dir}/nextflow_work \\
            ${ENSEMBL_ROOT_DIR}/ensembl-genomio/pipelines/nextflow/workflows/dumper_pipeline/main.nf \\
            -profile lsf \\
            \$(\${CMD} details script) \\
            --dbname_re '^drosophila_melanogaster_\\w+_57_.*\$' \\
            --output_dir \${data_dir}/dumper_output
        popd

        Mandatory arguments:
        --host, --port, --user           Connection parameters to the SQL servers we getting core db(s) from

        Optional arguments:
        --password                     Password part of the connection parameters
        --prefix                       Core dabase(s) name prefixes
        --dbname_re                    Regexp to match core db name(s) against
        --brc_mode	               Override Ensembl 'species' and 'division' with the corresponding BRC4 ones ('organism_abbrev' and 'component')
        --output_dir                   Name of Output directory to gather prepared outfiles. Default -> 'Output_GenomePrepare'.
        --help                         This usage statement.
        """
}

// Check mandatory parameters
if (params.help) {
    helpMessage()
    exit 0
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
        "prefix": "",
        "dbname_re": ""
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

if (params.host && params.port && params.user && params.output_dir) {
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
    DUMP_SQL(server, dbs, filter_map, params.output_dir)
    DUMP_METADATA(server, dbs, filter_map, params.output_dir)
}
