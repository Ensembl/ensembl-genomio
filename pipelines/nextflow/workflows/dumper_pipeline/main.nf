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

// Predefine the files that can be dumped and the number of files for each of them
default_selection_map = [
    'agp': 1,
    'annotation': 1,
    'seq_regions': 1,
    'events': 1,
    'genome_metadata': 1,
    'gff3': 1,
    'stats': 3,
    'fasta_pep': 1,
    'fasta_dna': 1,
    'seq_attrib': 1,
]
default_selection = default_selection_map.keySet() as ArrayList
params.db_list = ''
params.server_url = ''
params.host = ''

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
if (params.help) {
     log.info paramsHelp("nextflow run dumper_pipeline/main.nf --dump_sql --dump_all_files --host 'HOST' --port 'PORT' --user 'USER' --dbname_re 'DB_REGEX' --output_dir 'OUTPUT_DIR'")
    exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

if (params.brc_mode) {
    params.brc_mode = params.brc_mode as Integer
}

// Prepare the server params
include { extract_url_args; generate_url } from '../../modules/utils/utils.nf'
def create_server(params) {
    if (params.server_url) {
        server = extract_url_args(params.server_url)
        server.url = params.server_url
    } else if (params.host) {
        server = [
            "host": params.host,
            "port": params.port,
            "user": params.user,
        ]
        if (params.password) {
            server["password"] = params.password
        }
        server.url = generate_url("mysql", params.host, params.port, params.user, params.password)
    } else {
        log.info paramsHelp(
            "Server URL (--server_url) or parameters needed (--host, --port, --user, --password)"
        )
        exit 0
    }
    print("Using server: " + server)
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

server = create_server(params)
filter_map = create_filter_map(params)

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
    // Prepare db_factory filter
    basic_filter = Channel.of(filter_map)
    // Add db_list
    filter_db = Channel.of()
    if (params.db_list) {
        // Get the list of databases from the file
        db_list = Channel.fromPath(params.db_list).splitCsv().map{ row -> row[0] }.collect()
        // Add that list to the filter_db map
        filter_db = basic_filter.merge(db_list) { filters, db_list -> filters + ["db_list": db_list] }
    } else {
        filter_db = basic_filter
    }

    dbs = DB_FACTORY(server, filter_db)
        .map(it -> read_json(it))
        .flatten()
        .map(it -> it + ["id": it["production_name"]])
    
    if (dump_sql) {
        DUMP_SQL(dbs)
    }
    DUMP_FILES(dbs, dump_selection, dump_number)
}
