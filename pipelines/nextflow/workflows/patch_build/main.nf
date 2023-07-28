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

// Print usage
def helpMessage() {
  log.info '''
        Mandatory arguments:
        --host
        --port
        --user
        --pass                         Connection parameters to the SQL servers we getting core db(s) from
        --database                     Database name
        --events_file                  Annotation event file from gene_diff

        Optional arguments:
        --output_dir                   Where the process files can be stored
        --help                         This usage statement
        '''
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
        "password": params.pass,
        "database": params.database
    ]
    return server
}

if (params.host && params.port && params.user && params.pass) {
    server = create_server(params)
} else {
    exit 1, "Missing server parameters"
}

process extract_gene_lists {
    label 'local'

    input:
        path(events, stageAs: "events.tab")

    output:
        path("new_genes.tab")
        path("changed_genes.tab")

    script:
    def changed_genes = "changed_genes.tab"
    def new_genes = "new_genes.tab"
    """
    grep -v "//" $events | grep -v "=" | grep -v "~" | sed s'/[+><]/\\t/' | cut -f3 | sed 's/:/\\n/g' | sort -u > $new_genes
    grep "=" $events | cut -f2 | sed -r 's/=[+!-]?/\t/' > $changed_genes
    """
}

process publish {
    label "local"
    publishDir "$output_dir", mode: "copy"

    input:
    path(files)
    val(output_dir)

    output:
    path(files, includeInputs: true)
    
    script:
    """
    echo "$output_dir"
    """
}

// Run main workflow
workflow {
    events = Channel.fromPath(params.events_file, checkIfExists: true)

    // Extract the genes lists from the annotation event file
    (new_genes, changed_genes) = extract_gene_lists(events)

    // Temporary: get all generated files in one folder
    all_files = new_genes.concat(changed_genes)
    publish(all_files, params.output_dir)
}
