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
        --mock_osid                    Set to 1 if you want to run a fake id generator instead of OSID (for testing)
        or
        --osid_url
        --osid_user
        --osid_pass                    
        --osid_species                 Connection parameters to the OSID server

        Release parameters:
        --release_name                 Short release name (e.g. "66")
        --release_date                 Date of the release in ISO format (e.g. "2023-07-01")

        Optional arguments:
        --output_dir                   Where the logs and process files can be stored
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
        "database": params.database,
        "species": params.species
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

include { PATCH_BUILD_PROCESS } from '../../subworkflows/patch_build_post_process/main.nf'
include { PUBLISH_FILES } from '../../modules/patch_build/publish_files.nf'
workflow {
    events = Channel.fromPath(params.events_file, checkIfExists: true)
    deleted = Channel.fromPath(params.deletes_file, checkIfExists: true)
    old_registry = Channel.fromPath(params.old_registry, checkIfExists: true)
    new_registry = Channel.fromPath(params.new_registry, checkIfExists: true)
    server = create_server()
    osid = create_osid_params()
    release = create_release()

    logs = PATCH_BUILD_PROCESS(events, deleted, old_registry, new_registry, server, osid, release)
    if (params.output_dir) {
        PUBLISH_FILES(logs, params.output_dir)
    }
}
