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

// Patch build post-processing pipeline

// Import modules/subworkflows
include { CHECK_OSID } from '../../modules/patch_build/check_osid.nf'
include { EXTRACT_GENE_LISTS } from '../../modules/patch_build/extract_gene_lists.nf'
include { TRANSFER_IDS } from '../../modules/patch_build/transfer_ids.nf'
include { TRANSFER_METADATA } from '../../modules/patch_build/transfer_metadata.nf'
include { FORMAT_EVENTS } from '../../modules/patch_build/format_events.nf'
include { LOAD_EVENTS } from '../../modules/patch_build/load_events.nf'
include { ALLOCATE_IDS as ALLOCATE_GENE_IDS } from '../../modules/patch_build/allocate_ids.nf'
include { ALLOCATE_IDS as ALLOCATE_TRANSCRIPT_IDS } from '../../modules/patch_build/allocate_ids.nf'
include { FINALIZE_VERSIONS } from '../../modules/patch_build/finalize_versions.nf'
include { CHECK_PATCH } from '../../modules/patch_build/check_patch.nf'
include { CANONICAL_TRANSCRIPTS } from '../../modules/patch_build/canonical_transcripts.nf'

workflow PATCH_BUILD_PROCESS {

    take:
        events
        deleted
        old_registry
        new_registry
        server
        osid_params
        release
    
    emit:
        logs
    
    main:
        // Check OSID is available and has the expected species
        osid = CHECK_OSID(osid_params)

        // Extract the genes lists from the annotation event file
        (new_genes, changed_genes) = EXTRACT_GENE_LISTS(events)

        // Transfer the genes from the old db to the new db
        // Requires 2 registries and the species name
        transfer_ids_output = TRANSFER_IDS(changed_genes, old_registry, new_registry, server.species)
        new_transcripts = transfer_ids_output.new_transcripts
        mapped_transcripts = transfer_ids_output.mapped_transcripts

        // Transfer metadata (can be done any time after the ids are transfered?)
        transfer_log = TRANSFER_METADATA(mapped_transcripts, old_registry, new_registry, server.species)

        // Allocate ids for both the new_genes and the changed_genes new transcripts
        new_genes_map = ALLOCATE_GENE_IDS(new_registry, server.species, osid, new_genes, "gene")
        new_transcripts_map = ALLOCATE_TRANSCRIPT_IDS(new_registry, server.species, osid, new_transcripts, "transcript")

        // Finalize the versions
        waited_files = new_genes_map.concat(new_transcripts_map).last()
        finalized = FINALIZE_VERSIONS(server, waited_files)

        // Check stable ids
        patch_errors=CHECK_PATCH(server, finalized)

        // Regenerate the canonical transcripts
        CANONICAL_TRANSCRIPTS(server, waited_files)

        // Format the annotation events file into a compatible event file
        events_file = FORMAT_EVENTS(events, deleted, new_genes_map, release)
        LOAD_EVENTS(server, events_file)

        // Get out all generated files if we want to review them
        logs = changed_genes
            .concat(new_genes)
            .concat(new_genes_map)
            .concat(new_transcripts_map)
            .concat(new_transcripts)
            .concat(mapped_transcripts)
            .concat(events_file)
            .concat(transfer_log)
            .concat(patch_errors)
}
