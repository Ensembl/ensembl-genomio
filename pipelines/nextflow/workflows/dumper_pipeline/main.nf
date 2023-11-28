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
params.prefix = ''
params.dbname_re = ''
params.output_dir = './dumper_output'
params.password = ''

// Predefine the files that can be dumped and the number of files for each of them
params.dump_sql = false
params.dump_all_files = false
params.dump_selection = ''
default_selection_map = [
    'seq_regions': 1,
    'events': 1,
    'genome_metadata': 1,
    'stats': 3,
    'fasta_pep': 1,
    'fasta_dna': 1
]
default_selection = default_selection_map.keySet() as ArrayList

// Print usage
def helpMessage() {
  log.info """
        Mandatory arguments:
        --host, --port, --user         Connection parameters to the SQL servers we getting core db(s) from
        
        Dump SQL:
        --dump_sql                     Dump the whole database in SQL in a coredb subfolder

        Dump files:
        --dump_all_files               Dump all files (SQL excepted) in a metadata subfolder
        or
        --dump_selection               Dump files from a comma-separated list (among ${default_selection})

        Optional arguments:
        --password                     Password part of the connection parameters
        --prefix                       Core dabase(s) name prefixes
        --dbname_re                    Regexp to match core db name(s) against
        --brc_mode                     Override Ensembl 'species' and 'division' with the corresponding BRC ones ('organism_abbrev' and 'component')
        --output_dir                   Name of Output directory to gather prepared outfiles. (default: ${params.output_dir})
        --cache_dir                    Directory where some files are cached (e.g. NCBI stats files)
        --help                         This usage statement.

        Usage:
        The typical command for running the 'Dumper' pipeline is as follows:

        nextflow run \\
            -w \${data_dir}/nextflow_work \\
            ensembl-genomio/pipelines/nextflow/workflows/dumper_pipeline/main.nf \\
            --dump_sql --dump_files \\
            -profile lsf \\
            --host <DB_HOST> --port <DB_PORT> --user <DB_USER>
            --dbname_re '^drosophila_melanogaster_\\w+_57_.*\$' \\
            --output_dir \${data_dir}/dumper_output

        """
}

// Check mandatory parameters
if (params.help) {
    helpMessage()
    exit 0
}

if (params.brc_mode) {
    params.brc_mode = params.brc_mode as Integer
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

// Select the files to dump
dump_sql = false
dump_selection = []
if (params.dump_sql) {
    dump_sql = true
}
if (params.dump_all_files) {
    dump_selection = default_selection
} else if (params.dump_selection) {
    dump_selection = params.dump_selection.split(/,/).collect().unique()
    for (item in dump_selection) {
        if (!default_selection.contains(item)) {
            acceptable = default_selection.join(", ")
            exit 1, "Selection item unknown: " + item + " (accepted: " + acceptable + ")"
        }
    }
}

// Compute the number of files to dump
dump_number = 0
for (dump_key in dump_selection) {
    dump_number += default_selection_map[dump_key]
}
if (!dump_sql and dump_number == 0) {
    exit 1, "No dump option selected"
}
if (dump_sql) {
    print("Will dump databases to SQL in $params.output_dir")
}
if (dump_number) {
    print("Will dump $dump_number files per genome for $dump_selection in $params.output_dir")
}

include { DUMP_SQL } from '../../subworkflows/dump_sql/main.nf'
include { DUMP_FILES } from '../../subworkflows/dump_files/main.nf'
include { DB_FACTORY } from '../../modules/database/db_factory.nf'
include { read_json } from '../../modules/utils/utils.nf'

// Run main workflow
workflow {
    dbs = DB_FACTORY(server, filter_map)
        .map(it -> read_json(it))
        .flatten()
    
    if (dump_sql) {
        DUMP_SQL(dbs)
    }
    DUMP_FILES(dbs, dump_selection, dump_number)
}
