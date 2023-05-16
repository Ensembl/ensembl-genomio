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

// Import modules/subworkflows
include { download_genbank; extract_from_gb} from '../modules/download_genbank.nf'
include { process_gff3; gff3_validation } from '../modules/process_gff3.nf'
include { CHECK_JSON_SCHEMA } from '../modules/check_json_schema.nf'
include { JSON_SCHEMA_FACTORY } from '../modules/json_schema_factory.nf'
include { COLLECT_FILES; MANIFEST; PUBLISH_DIR } from '../modules/collect_files.nf'
include { CHECK_INTEGRITY } from '../modules/integrity.nf'
include { MANIFEST_STATS } from '../modules/manifest_stats.nf'

workflow additional_seq_prepare {
    take:
        prefix
        accession
        production_name
        brc_mode
    emit:
        manifeststats_ch 
    main:
        gb_file_ch = download_genbank(accession)
        extract_from_gb(prefix, production_name, gb_file_ch)
        annotation = process_gff3(extract_from_gb.out.gene_gff, extract_from_gb.out.genome, accession)
        gff3_validation(process_gff3.out.gene_models)
        json_files = extract_from_gb.out.genome.concat(extract_from_gb.out.seq_regions, process_gff3.out.functional_annotation)
        json_files_checked = CHECK_JSON_SCHEMA(json_files, annotation.out.gca)
                        .map{verified_json, gca -> verified_json}
        all_files = json_files_checked
                        .concat(extract_from_gb.out.dna_fasta, extract_from_gb.out.pep_fasta)
                        .map{it -> [accession, it]}
                        .groupTuple()
        collect_dir = COLLECT_FILES(all_files)
        manifested_dir = MANIFEST(collect_dir, accession)
        manifest_checked = CHECK_INTEGRITY(manifested_dir, brc_mode)
        PUBLISH_DIR(manifest_checked, accession)
        manifeststats_ch = MANIFEST_STATS(manifested_dir, accession, 'datasets')
}

workflow {
    manifest_stats = additional_seq_prepare(params.prefix, params.accession, params.production_name, params.brc_mode)
    manifest_stats.view()
}
