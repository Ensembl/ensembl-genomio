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
params.mock_osid = false

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

        Registries parameters:
        --old_registry                 Registry where the old version of the db resides (transfer from)
        --new_registry                 Registry where the new version of the db resides (to modify) 
        --species                      Production name of the species, same in both registries
        
        OSID parameters:
        --osid_url
        --osid_user
        --osid_pass                    
        --osid_species                 Connection parameters to the OSID server
        --mock_osid                    Set to 1 if you want to run a fake id generator instead of OSID (for testing)

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

if (!params.host or !params.port or !params.user or !params.pass) {
    exit 1, "Missing server parameters"
}
if (!params.mock_osid) {
    if (!params.osid_url or !params.osid_user or !params.osid_pass or !params.osid_species) {
        exit 1, "Missing OSID parameters"
    }
}
if (!params.old_registry or !params.new_registry or !params.species) {
    exit 1, "Missing registries parameters"
}

def create_server() {
    server = [
        "host": params.host,
        "port": params.port,
        "user": params.user,
        "password": params.pass,
        "database": params.database
    ]
    return server
}

def create_osid_params() {
    if (params.mock_osid) {
        osid = ["mock": true]
    } else {
        osid = [
            "url": params.osid_url,
            "user": params.osid_user,
            "pass": params.osid_pass,
            "species": params.osid_species
        ]
    }
    return osid
}

def create_osid_server_params() {

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

process transfer_ids {
    label 'local'

    input:
        path(changed_genes)
        path(old_registry)
        path(new_registry)
        val(species)

    output:
        path("new_transcripts.txt")

    script:
    def new_transcripts = "new_transcripts.txt"
    """
    touch $new_transcripts
    """
}

process transfer_metadata {
    label 'local'

    input:
        path(old_registry)
        path(new_registry)
        val(species)

    script:
    """
    """
}

process allocate_ids {
    label 'local'

    input:
        path(new_registry)
        val(species)
        val(osid_params)
        path(new_genes)
        val(mode)
    
    output:
        path("output_map.txt")

    script:
    def output_map = "output_map.txt"
    """
    touch $output_map
    """
}

process format_events {
    label 'local'

    input:
        path(events, stageAs: "annotation_events.txt")
        path(new_genes)
    
    output:
        path("events.tab")

    script:
    def final_events = "events.tab"
    """
    touch $final_events
    """
}

process load_events {
    label 'local'

    input:
        val(server)
        path(events)

    script:
    """
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


workflow PATCH_BUILD_PROCESS {

    take:
        events
        old_registry
        new_registry
        server
        osid
    
    main:
        // Extract the genes lists from the annotation event file
        (new_genes, changed_genes) = extract_gene_lists(events)

        // Transfer the genes from the old db to the new db
        // Requires 2 registries and the species name
        new_transcripts = transfer_ids(changed_genes, old_registry, new_registry, params.species)

        // Transfer metadata (can be done any time after the ids are tranfered?)
        transfer_metadata(old_registry, new_registry, params.species)

        // Allocate ids for both the new_genes and the changed_genes new transcripts
        new_genes_map = allocate_ids(new_registry, params.species, osid, new_genes, "gene")
        // new_transcripts_map = allocate_transcript_ids(new_registry, species, osid_params, new_transcripts)

        // Format the annotation events file into a compatible event file
        events_file = format_events(events, new_genes_map)
        load_events(server, events_file)

        // Temporary: get all generated files in one folder
        all_files = new_genes
            .concat(changed_genes)
            .concat(new_transcripts)
            .concat(new_genes_map)
            .concat(events_file)
        publish(all_files, params.output_dir)
}

// Run main workflow
workflow {
    events = Channel.fromPath(params.events_file, checkIfExists: true)
    old_registry = Channel.fromPath(params.old_registry, checkIfExists: true)
    new_registry = Channel.fromPath(params.new_registry, checkIfExists: true)
    server = create_server()
    osid = create_osid_params()

    PATCH_BUILD_PROCESS(events, old_registry, new_registry, server, osid)
}
