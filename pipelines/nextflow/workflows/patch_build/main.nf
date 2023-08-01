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

        Input files parameters:
        --events_file                  Annotation event file from gene_diff
        --deletes_file                 Deleted genes file

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

        Release parameters:
        --release_name                 Short release name (e.g. "66")
        --release_date                 Date of the release in ISO format (e.g. "2023-07-01")

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

if (!params.scripts_dir) {
    exit 1, "Missing scripts_dir"
}

if (!params.events_file or !params.deletes_file) {
    exit 1, "Missing input files parameters"
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
if (!params.release_name or !params.release_date) {
    exit 1, "Missing release parameters"
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

def create_release() {
    release = [
        "name": params.release_name,
        "date": params.release_date,
    ]
    return release
}

include { EXTRACT_GENE_LISTS } from '../../modules/patch_build/extract_gene_lists.nf'
include { TRANSFER_IDS } from '../../modules/patch_build/transfer_ids.nf'
include { TRANSFER_METADATA } from '../../modules/patch_build/transfer_metadata.nf'
include { FORMAT_EVENTS } from '../../modules/patch_build/format_events.nf'
include { LOAD_EVENTS } from '../../modules/patch_build/load_events.nf'
include { PUBLISH_FILES } from '../../modules/patch_build/publish_files.nf'
include { ALLOCATE_IDS as ALLOCATE_GENE_IDS } from '../../modules/patch_build/allocate_ids.nf'
include { ALLOCATE_IDS as ALLOCATE_TRANSCRIPT_IDS } from '../../modules/patch_build/allocate_ids.nf'

workflow PATCH_BUILD_PROCESS {

    take:
        events
        deleted
        old_registry
        new_registry
        server
        osid
        release
    
    main:
        // Extract the genes lists from the annotation event file
        (new_genes, changed_genes) = EXTRACT_GENE_LISTS(events)

        // Transfer the genes from the old db to the new db
        // Requires 2 registries and the species name
        new_transcripts = TRANSFER_IDS(changed_genes, old_registry, new_registry, params.species)

        // Transfer metadata (can be done any time after the ids are transfered?)
        transfer_log = TRANSFER_METADATA(new_transcripts, old_registry, new_registry, params.species)

        // Allocate ids for both the new_genes and the changed_genes new transcripts
        new_genes_map = ALLOCATE_GENE_IDS(new_registry, params.species, osid, new_genes, "gene")
        new_transcripts_map = ALLOCATE_TRANSCRIPT_IDS(new_registry, params.species, osid, new_transcripts, "transcript")

        // Format the annotation events file into a compatible event file
        events_file = FORMAT_EVENTS(events, deleted, new_genes_map, release)
        LOAD_EVENTS(server, events_file)

        // Temporary: get all generated files in one folder
        all_files = changed_genes
            .concat(new_genes)
            .concat(new_genes_map)
            .concat(new_transcripts)
            .concat(events_file)
            .concat(transfer_log)
        PUBLISH_FILES(all_files, params.output_dir)
}

// Run main workflow
workflow {
    events = Channel.fromPath(params.events_file, checkIfExists: true)
    deleted = Channel.fromPath(params.deletes_file, checkIfExists: true)
    old_registry = Channel.fromPath(params.old_registry, checkIfExists: true)
    new_registry = Channel.fromPath(params.new_registry, checkIfExists: true)
    server = create_server()
    osid = create_osid_params()
    release = create_release()

    PATCH_BUILD_PROCESS(events, deleted, old_registry, new_registry, server, osid, release)
}
