#!env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import eHive
import gzip
import io
import json
import os
import re
import subprocess as sp
import sys

from collections import defaultdict
from os.path import dirname, join as pj


class load_sequence_data(eHive.BaseRunnable):
    """
    loading sequence data, seq region names, atrributes and synonyms

    eHive module to load sequnce data, seq region names, atrributes and synonyms from FASTAs AGPs and seq_region.json.
    Various ensembl-analysis perl scripts are used to create coord_systems, load sequences and set attributes.
    SQL commands through out the code to be replaces with the proper python API at some point.
    """

    def param_defaults(self):
        """
        default parameter/options values
        """
        return {
            # relative order of the coord_systems types to infer their rank from
            "cs_order": "ensembl_internal,chunk,contig,supercontig,non_ref_scaffold,scaffold,primary_assembly,superscaffold,linkage_group,chromosome",
            "artificial_cs": "ensembl_internal,chunk",  # coord_systems to ignore when adding synonyms for sequence_level cs
            "IUPAC": "RYKMSWBDHV",  # symbols to be replaced with N in the DNA sequences (ensembl core(107) doesn't support the whole IUPAC alphabet for DNA)
            # unversion scaffold, remove ".\d$" from seq_region.names if there's a need
            "unversion_scaffolds": 0,
            "versioned_sr_syn_src": "INSDC",  # INSDC(50710) # if unversioning non-sequence level cs, store original name (with version) as this synonym
            "sr_syn_src": "BRC4_Community_Symbol",  # BRC4_Community_Symbol(211) # if unversioning sequence-level cs, store original name (with version) as this synonym
            # nullify coord_system version for the given coord_system name
            "nullify_cs_version_from": "contig",
            # default coord_system name for single-level (no AGPs assemblies)
            "noagp_cs_name_default": "primary_assembly",
            # file for additional mapping of the synonym sources to the external_db (ensembl), as used by "get_external_db_mapping" function below
            "external_db_map": None,
            # set "coord_system_tag" seq_region attribute to this value if there's a corresponding "chromosome_display_order" list in genome.json metadata
            #   if None, only "chromosome" coord system is processed (if present)
            #   (see add_chr_karyotype_rank definition below )
            #   if seq_region already has `coord_system_tagi` attribute, it value updated only if "force_update_coord_system_tag" module param is True (see below)
            "cs_tag_for_ordered": None,
            # Force updating of the "coord_system_tags" attribute for seq_regions from `cs_tag_for_ordered` (see above)
            "force_update_coord_system_tag": False,
            # BRC4 compatibility mode; if on, "(EBI|BRC4)_seq_region_name" seq_region_attributes are added.
            #   Blocked by the "swap_gcf_gca" option. In this case insertion should be done on later pipeline stage after seq_region name swapping.
            "brc4_mode": True,
            # Whether to use RefSeq names as additional seq_region synonyms (if available) or not (see add_sr_synonyms definition below)
            #  Does not swap anything actually, just loads synonyms to be used by a later "swapping" stage
            #  Disables BRC4 compatibilty mode (see the "brc4_mode" option comment).
            "swap_gcf_gca": False,
            # list of coord systems used in "-ignore_coord_system" options of the "ensembl-analysis/scripts/assembly_loading/set_toplevel.pl" script
            #   part of the loading process
            "not_toplevel_cs": [],  # i.e. "contig", "non_ref_scaffold"
            # explicit list of seq_region properties (keys) to load as seq_region_attrib s (values) (see add_sr_attribs definition below)
            #   if a dict's used as a value, treat its keys as "json_path" (/ as delim) map, i.e.
            #       { "added_sequence" : { "assembly_provider" : { "name" : ... } } } -> "added_sequence/assembly_provider/name"
            #   only flattable properties can be used, no arrays
            #   arrays should be processed separately (see `add_sr_synonyms` or `add_karyotype_bands` definitions)
            # see src/python/ensembl/io/genomio/data/schemas/seq_region.json
            "sr_attrib_types": {
                "circular": "circular_seq",
                "codon_table": "codon_table",
                "location": "sequence_location",
                "non_ref": "non_ref",
                "coord_system_level": "coord_system_tag",
                "added_sequence": {
                    # json_path to attrib_type_code(str) mapping
                    "added_sequence/accession": "added_seq_accession",
                    "added_sequence/assembly_provider/name": "added_seq_asm_pr_nam",
                    "added_sequence/assembly_provider/url": "added_seq_asm_pr_url",
                    "added_sequence/annotation_provider/name": "added_seq_ann_pr_nam",
                    "added_sequence/annotation_provider/url": "added_seq_ann_pr_url",
                },  # added_sequence
            },
            # loading additional sequences to the already exsisting core db
            "load_additional_sequences": 0,
            # size of the sequence data chunk, if 0 (default), no chunking is performed
            "sequence_data_chunk": 0,
            #   min size of the sequence chunk, no chunking is done if 'sequence_data_chunk' < 'sequence_data_chunk_min'
            "sequence_data_chunk_min_len": 50_000,
            # coord system name for chunks
            "chunk_cs_name": "ensembl_internal",
        }

    def run(self):
        """
        Entry point for the Ehive module. All processing is done here in this case.
        """
        # params
        work_dir = self.param_required("work_dir")

        # initial sequence loading, using ensembl-analysis scripts
        self.initial_sequence_loading(work_dir)

        # load data from the corresponding core db tables
        external_db_map = self.load_map_from_core_db(
            "external_db", ["db_name", "external_db_id"], work_dir
        )  # for external_db
        attrib_type_map = self.load_map_from_core_db(
            "attrib_type", ["code", "attrib_type_id"], work_dir
        )  # for attrib_type
        seq_region_map = self.load_map_from_core_db(
            "seq_region", ["name", "seq_region_id"], work_dir
        )  # for seq_region

        # update synonyms and seq_region_attribs
        unversion = self.param("unversion_scaffolds")
        is_primary_assembly = self.from_param("manifest_data", "agp", not_throw=True) is None
        seq_region_file = self.from_param("manifest_data", "seq_region", not_throw=True)

        #   add seq_region synonyms
        self.add_sr_synonyms(
            seq_region_file,
            seq_region_map,
            external_db_map,
            self.pjc(work_dir, "seq_region_syns"),
            unversion=unversion,
        )

        #   add seq_region attributes
        self.add_sr_attribs(
            seq_region_file,
            seq_region_map,
            attrib_type_map,
            self.pjc(work_dir, "seq_region_attr"),
            unversion=unversion,
        )

        #   add seq_region EBI and BRC4 name attributes in the "BRC4 mode"
        #     special case of attributes adding with default values derived from seq_region names
        #     do not add if preparing to swap RefSeq and GeneBank ids; in this case attributes to be added at a later stage in pipeline
        #     (easier to insert then to update)
        if self.param("brc4_mode") and not self.param("swap_gcf_gca"):
            self.add_sr_ebi_brc4_names(
                seq_region_file,
                seq_region_map,
                attrib_type_map,
                self.pjc(work_dir, "seq_region_ebi_brc4_name"),
                unversion=unversion,
            )

        # add karyotype related data
        self.add_karyotype_data(
            seq_region_file,
            seq_region_map,
            attrib_type_map,
            self.pjc(work_dir, "karyotype"),
            unversion=unversion,
        )

    def initial_sequence_loading(self, work_dir: str):
        """
        initial preparation and loading of AGPs and fasta data.

        initial preparation and loading of AGPs and fasta data using ensembl-analysis perl scripts
        """
        # preprocess FASTA with sequences
        #   rename IUPAC to N symbols using sed
        fasta_clean = self.from_param("manifest_data", "fasta_dna")

        # start coord system ranking and agps processing
        agps = self.from_param("manifest_data", "agp", not_throw=True)

        # get the deafult coord_system order
        #   use noagp_cs_name_default for "noagp" assemblies
        cs_order = self.coord_sys_order(self.param("cs_order"))
        noagps_cs = self.param("noagp_cs_name_default")

        # remove gaps and lower_level mappings if the are coveres by higher level ones
        #   i.e.: remove 'contigN to chromosomeZ', if 'contigN to scaffoldM' and 'scaffoldM to chromosomeZ' are in place
        #   returns None if no agps provided
        agps_pruned_dir = self.pjc(work_dir, "agps_pruned")
        agps_pruned = self.prune_agps(agps, cs_order, agps_pruned_dir, self.param_bool("prune_agp"))

        # order
        # rank cs_names, met in agps.keys ("-" separated, i.e. "scaffold-contig") based on cs_order
        cs_rank = self.used_cs_ranks(agps_pruned, cs_order, noagps_cs)

        # chunk sequence data if needed
        #   no chunking if chunk_size < 50k
        chunk_size = int(self.param("sequence_data_chunk"))
        chunk_cs_name = self.param("chunk_cs_name")
        fasta_clean, cs_rank, agps_pruned = self.chunk_contigs(
            fasta_clean,
            cs_rank,
            agps_pruned,
            pj(work_dir, "chunking"),
            chunk_size=chunk_size,
            chunks_cs_name=chunk_cs_name,
        )

        # empty agps_pruned ignored
        self.load_seq_data(fasta_clean, agps_pruned, cs_rank, self.pjc(work_dir, "load"))

        # mark all the "contig"s or noagp_cs as being sourced from ENA
        if not self.param_bool("no_contig_ena_attrib"):
            # NB using original "agps" parameter (with no chunking data added)
            agps_raw = self.from_param("manifest_data", "agp", not_throw=True)
            if agps_raw is None:
                self.add_contig_ena_attrib(self.pjc(work_dir, "load", "set_ena"), cs_name=noagps_cs)
            else:
                self.add_contig_ena_attrib(self.pjc(work_dir, "load", "set_ena"))

        # unversion scaffold, remove ".\d$" from names if there's a need
        if self.param_bool("unversion_scaffolds"):
            self.unversion_scaffolds(cs_rank, self.pjc(work_dir, "unversion_scaffolds"))

        # add assembly mappings between various cs to meta table for the mapper to work properly
        cs_pairs = agps_pruned and agps_pruned.keys() or None
        self.add_asm_mappings(cs_pairs, self.pjc(work_dir, "asm_mappings"))

        # set toplevel seq_region attribute
        self.set_toplevel(self.pjc(work_dir, "set_toplevel"), self.param("not_toplevel_cs"))

        # nullify contig version and update mappings strings accordingly; ignore for "load_additional_sequences" mode
        if not self.param_bool("load_additional_sequences"):
            self.nullify_ctg_cs_version(cs_order, self.pjc(work_dir, "asm_mapping", "nullify_cs_versions"))

    def add_sr_synonyms(
        self,
        seq_region_file: str,
        seq_region_map: dict,
        external_db_map: dict,
        work_dir: str,
        unversion: bool = False,
        unversionable_sources_set: frozenset = frozenset(["INSDC", "RefSeq"]),
    ):
        """
        Add seq_region_synonym from the seq_region_file meta data file.

        Add seq_region_synonym from the src/python/ensembl/io/genomio/data/schemas/seq_region.json compatible meta data file.
        Merge with the already existing ones in the db.

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible
          * the unversioned synonyms from the unversionable_sources_set will be added as well as the original ones

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file:
            return

        # get seq_region ids, names, syns from db
        synonyms_trios_db = self.load_seq_region_synonyms_trios_from_core_db(
            self.pjc(work_dir, "syns_from_core")
        )

        # form set of synonyms already present in db
        synonyms_in_db = frozenset([trio[2] for trio in synonyms_trios_db if trio[2] != "NULL"])

        # subset of sources to use the allowed the unversion synonyms
        unversionable_sources = unversion and unversionable_sources_set or frozenset()

        # technical / optimization. get external_db_id for "ensembl_internal_synonym"
        ensembl_internal_synonym_ext_db_id = self.id_from_map_or_die(
            "ensembl_internal_synonym", external_db_map, "external_db_map"
        )

        # get dict for additional mapping for sources to external_db names if there's one specified by "external_db_map" module param
        #   not to be confused with the external_db_map function parameter above
        additional_sources_mapping = self.get_external_db_mapping()

        # load synonyms from the json file
        synonyms_from_json = (
            []
        )  # [ (seq_region_id, synonym, external_db_id)... ] list of trios for inserting into db
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in filter(lambda sr: sr.get("synonyms", False), seq_regions):
                # iterate through all seq_regions having "synonyms"
                seq_region_name, seq_region_id, _ = self.name_and_id_from_seq_region_item(
                    seq_region, seq_region_map, try_unversion=unversion
                )

                # fill synonyms_from_json list of trios
                for synonym_item in list(seq_region["synonyms"]):
                    synonym_name = synonym_item["name"]
                    source = synonym_item["source"]
                    unversioned_name = ""

                    # check if there's any addtional mapping for the source, remap if so
                    if additional_sources_mapping:
                        source = additional_sources_mapping.get(
                            source, source
                        )  # use the same name if no matches in additional_sources_mapping dict

                    # try to get unversioned name if applicable
                    if source in unversionable_sources:
                        unversioned_name = re.sub(r"\.\d+$", "", synonym_name)

                    # put trios if names are not already seen in db
                    if synonym_name not in synonyms_in_db:
                        external_db_id = self.id_from_map_or_die(source, external_db_map, "external_db_map")
                        synonyms_from_json.append(
                            (seq_region_id, self.quote_or_null(synonym_name), external_db_id)
                        )

                    #   put additional unversioned synonyms if there's a sane one
                    if (
                        unversioned_name
                        and unversioned_name != synonym_name
                        and unversioned_name not in synonyms_in_db
                    ):
                        synonyms_from_json.append(
                            (
                                seq_region_id,
                                self.quote_or_null(unversioned_name),
                                ensembl_internal_synonym_ext_db_id,
                            )
                        )

        # run insertion SQL
        self.insert_to_db(
            synonyms_from_json,
            "seq_region_synonym",
            ["seq_region_id", "synonym", "external_db_id"],
            self.pjc(work_dir, "new_seq_region_synonyms"),
            ignore=True,
        )

    def add_sr_attribs(
        self,
        seq_region_file: str,
        seq_region_map: dict,
        attrib_type_map: dict,
        work_dir,
        unversion: bool = False,
    ):
        """
        Add seq_region_attrib(s) from the seq_region_file meta data file.

        Explicit list is taken from "sr_attrib_types" module param.

        Add seq_region_attrib(s) from the src/python/ensembl/io/genomio/data/schemas/seq_region.json compatible meta data file.
        Explicit list is taken from "sr_attrib_types" module param.

        "sr_attrib_types" defines { json_property -> attrib_type.name } map. If the value is dict,
        its keys are treated as "/"-delimited "json_path" (i.e. "added_sequence/assembly_provider/name").
        No arrays can be processed. Only simple or "flattable" types.

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file:
            return

        # technical / optimization. get atttib_type_id(s)
        # create a smaller map with attrib_type_id(s) as values
        properties_to_use = (
            []
        )  # [frozen]set with the top-level "seq_region" properties, that should be processed
        path_attrib_id_map = dict()  # { "flatterned/json/paths" : attrib_id_map ))}
        # fill set and map
        for prop, attrib_type in self.param("sr_attrib_types").items():
            # adding high level properties to process
            properties_to_use.append(prop)
            # adding json paths (or properties themselves) to atrrib_type_id map
            if isinstance(attrib_type, dict):  # if using json paths (delimeterd with "/")
                for path, inner_attrib_type in attrib_type.items():
                    path_attrib_id_map[path] = self.id_from_map_or_die(
                        inner_attrib_type, attrib_type_map, "attrib_type_map"
                    )
            else:
                path_attrib_id_map[prop] = self.id_from_map_or_die(
                    attrib_type, attrib_type_map, "attrib_type_map"
                )
        # return if there's nothing to add
        if not properties_to_use:
            return

        properties_to_use = frozenset(properties_to_use)

        # load attributes from seq_region file
        attrib_trios = []  # [ (seq_region_id, attrib_id, value)... ] list of trios for inserting into db
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in seq_regions:
                # get seq_region_id (perhaps, by using unversioned name)
                seq_region_name, seq_region_id, unversioned_name = self.name_and_id_from_seq_region_item(
                    seq_region, seq_region_map, try_unversion=unversion
                )

                # iterate through properties
                for prop_name in properties_to_use:
                    if prop_name not in seq_region:
                        continue
                    # flattern path
                    path_attrib_id_values_list = self.flattern_seq_region_item(
                        seq_region, prop_name, path_attrib_id_map
                    )
                    # fill attrib_trios
                    for path, attrib_id, value in path_attrib_id_values_list:
                        attrib_trios.append((seq_region_id, attrib_id, self.quote_or_null(value)))

        # run insertion SQL
        self.insert_to_db(
            attrib_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            self.pjc(work_dir, "brc4_ebi_seq_region_synonyms"),
            ignore=True,
        )

    def flattern_seq_region_item(
        self, seq_region: dict, prop_name: str, path_attrib_id_map: dict, sep: str = "/"
    ) -> list:
        """
        Flattern seq_region[property] and store corresponding [ (json_path, attrib_id, value)... ] (as list of trios).

        Only works for simple properties or dicts with no arrays on the path. Basically, implemets tree traversal.
        Utility function used by the `add_sr_attribs` method
        """
        res = []
        # is there anything to do
        if prop_name not in seq_region:
            return res

        # set up
        value = seq_region[prop_name]
        paths_to_go = [
            (prop_name, value)
        ]  # storing path and the corresponding value, to prevent repetetive traversals
        # iterate
        while paths_to_go:
            (path, value) = paths_to_go.pop()  # get last item
            if isinstance(value, list):
                # perhaps, it's better to raise exception then to continue silently
                continue
            if isinstance(value, dict):
                # if value is a complex object, add its leaves
                for key, val in value.items():
                    paths_to_go.append((f"{path}{sep}{key}", val))
                continue
            # if value is simple
            attrib_id = path_attrib_id_map.get(path, None)
            if attrib_id:
                res.append((path, attrib_id, value))
        # return what ever we have
        return res

    def add_sr_ebi_brc4_names(
        self,
        seq_region_file: str,
        seq_region_map: dict,
        attrib_type_map: dict,
        work_dir: str,
        unversion: bool = False,
    ):
        """
        Add "(EBI|BRC4)_seq_region_name" seq_region_attrib(s) either from the seq_region_file meta data file, or from original seq_region names.

        Add "(EBI|BRC4)_seq_region_name" seq_region_synonym from the src/python/ensembl/io/genomio/data/schemas/seq_region.json compatible meta data file or from the original seq_region_names.
        A special case of attributes adding with default values derived from seq_region names.

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file:
            return

        # technical / optimization. get atttib_type_id(s) for "(EBI|BRC4)_seq_region_name"
        tagged_sr_name_attrib_id = {
            tag: self.id_from_map_or_die(f"{tag}_seq_region_name", attrib_type_map, "attrib_type_map")
            for tag in ["EBI", "BRC4"]
        }

        # load BRC4/EBI name from seq_region file
        brc4_ebi_name_attrib_trios = (
            []
        )  # [ (seq_region_id, attrib_id, value)... ] list of trios for inserting into db
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in seq_regions:
                # get seq_region_id (perhaps, by using unversioned name)
                seq_region_name, seq_region_id, unversioned_name = self.name_and_id_from_seq_region_item(
                    seq_region, seq_region_map, try_unversion=unversion
                )
                # append attribs to the brc4_ebi_name_attrib_trios list
                for tag in ["BRC4", "EBI"]:
                    attrib_name = f"{tag}_seq_region_name"
                    attrib_id = tagged_sr_name_attrib_id[tag]
                    value = seq_region.get(attrib_name, seq_region_name)
                    brc4_ebi_name_attrib_trios.append((seq_region_id, attrib_id, self.quote_or_null(value)))

        # run insertion SQL
        self.insert_to_db(
            brc4_ebi_name_attrib_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            self.pjc(work_dir, "brc4_ebi_seq_region_synonyms"),
            ignore=True,
        )

    def add_karyotype_data(
        self,
        seq_region_file: str,
        seq_region_map: dict,
        attrib_type_map: dict,
        work_dir: str,
        unversion: bool = False,
    ):
        """
        Adds various karyotypic data from seq_region file and assembly metadata (if present).

        Adds various karyotypic data from the src/python/ensembl/io/genomio/data/schemas/seq_region.json compatible meta data file and assembly metadata (if present).

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible
        """
        # add karyotyope bands data
        regions_with_karyotype_bands = self.add_karyotype_bands(
            seq_region_file,
            seq_region_map,
            attrib_type_map,
            self.pjc(work_dir, "karyotype_bands"),
            unversion=unversion,
        )

        # try to add karyotype ranks for regions listed in genome_data/assembly/chromosome_display_order metadata
        regions_with_ranks_from_assembly_metadata = self.add_karyotype_rank_based_on_assembly_metadata(
            seq_region_map,
            attrib_type_map,
            self.pjc(work_dir, "karyotype_ranks_from_meta"),
            unversion=unversion,
        )

        regions_with_ranks_from_chromosome_cs = []
        if not regions_with_ranks_from_assembly_metadata:
            # try to add karyotype_ranks for top-level regions from the "chromosome" coord_system
            regions_with_ranks_from_chromosome_cs = self.add_karyotype_rank_for_chromosomes(
                attrib_type_map, self.pjc(work_dir, "karyotype_ranks_for_chromosomes")
            )

        # make sure that regions with bands have karyotype_ranks
        self.add_karyotype_rank_from_bands_info(
            regions_with_karyotype_bands,
            regions_with_ranks_from_chromosome_cs + regions_with_ranks_from_assembly_metadata,
            attrib_type_map,
            self.pjc(work_dir, "karyotype_ranks_from_bands"),
        )

    def add_karyotype_bands(
        self,
        seq_region_file: str,
        seq_region_map: dict,
        attrib_type_map: dict,
        work_dir: str,
        unversion: bool = False,
        karyotype_bands_property="karyotype_bands",
    ) -> list:  # [ (seq_region_name, seq_region_id, unversioned_name) ]
        """
        Add karyotypic data from the seq_region metafile.

        Add karyotypic data from the src/python/ensembl/io/genomio/data/schemas/seq_region.json compatible meta data file.
        Returns list of [ (seq_region_name, seq_region_id, unversioned_name) ] trios for seq_regions having karyotype bands info.

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file:
            return

        # resulting list of seq regions with bands
        seq_regions_with_karyotype_bands = []  # [ ( seq_region_name, seq_region_id, unversioned_name )... ]

        # load BRC4/EBI name from seq_region file
        band_tuples = (
            []
        )  # [ (seq_region_id, seq_region_start, seq_region_end, band|"NULL", stain|"NULL")... ] list of tuples for inserting into db
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in filter(lambda sr: sr.get(karyotype_bands_property, False), seq_regions):
                # iterate through all seq_regions having non-empty "karyotype_bands"

                # get seq_region_id (perhaps, by using unversioned name)
                seq_region_name, seq_region_id, unversioned_name = self.name_and_id_from_seq_region_item(
                    seq_region, seq_region_map, try_unversion=unversion
                )

                # append trio to the resulting list
                seq_regions_with_karyotype_bands.append((seq_region_name, seq_region_id, unversioned_name))

                # append bands to the band_tuples list
                for band in seq_region[karyotype_bands_property]:
                    # print("BAND: " + str(band), file = sys.stderr)
                    # coords
                    seq_region_start = band["start"]
                    seq_region_end = band["end"]
                    # band_name and stain
                    band_name = band.get("name", None)
                    stain = band.get("stain", None)
                    # special cases for stain
                    structure = band.get("structure", None)
                    if structure == "telomere":
                        stain = "TEL"
                    elif structure == "centromere":
                        stain = "ACEN"

                    # append tuple
                    band_tuples.append(
                        (
                            seq_region_id,
                            seq_region_start,
                            seq_region_end,
                            self.quote_or_null(band_name),
                            self.quote_or_null(stain),
                        )
                    )

        # run insertion SQL
        self.insert_to_db(
            band_tuples,
            "karyotype",
            ["seq_region_id", "seq_region_start", "seq_region_end", "band", "stain"],
            self.pjc(work_dir, "karyotype_insertion"),
            ignore=True,
        )

        # return resulting list of regions with bands trios
        return seq_regions_with_karyotype_bands

    def add_karyotype_rank_based_on_assembly_metadata(
        self, seq_region_map: dict, attrib_type_map: dict, work_dir: str, unversion: bool = True
    ) -> list:  # [ (seq_region_name, seq_region_id, unversioned_name) ]
        """
        Add `karyotype_rank` attributes for seq region data from based on metadata from the "genome_data" module parameter.
        Add only to the seq_regions with ids listed in the array corresponding to 'genome_data/assembly/chromosome_display_order'.

        Set "coord_system_tag" attribute to the one listed in the "cs_tag_for_ordered" module param; or "chromosome" if param value is underfined.
        Force updating of the "coord_system_tags" if `force_update_coord_system_tag` module param is True.

        Returns list of [ (seq_region_name, seq_region_id, unversioned_name) ] trios for seq_regions with updated karyotype_ranks.

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """

        # resulting list of seq_region with karyotype_rank
        regions_with_ranks_from_assembly_metadata = (
            []
        )  # [ ( seq_region_name, seq_region_id, unversioned_name )... ]

        # get `chromosome_display_order` list from the assembly metadata
        assembly_metadata = self.from_param("genome_data", "assembly", not_throw=True) or dict()
        chromosome_display_order_list = assembly_metadata.get("chromosome_display_order", [])

        # technical / optimization. get external_db_id for "karyotype_rank" and "coord_system_tag"
        karyotype_rank_attrib_id = self.id_from_map_or_die(
            "karyotype_rank", attrib_type_map, "attrib_type_map"
        )
        coord_system_tag_attrib_id = self.id_from_map_or_die(
            "coord_system_tag", attrib_type_map, "attrib_type_map"
        )
        coord_system_tag = self.param("cs_tag_for_ordered") or "chromosome"
        force_update_coord_system_tag = self.param("force_update_coord_system_tag") or False

        # set/update proper attributes for `chromosome_display_order_list` regions
        rank_insertions_trios = []
        coord_system_tag_attrib_insertion_trios = []
        coord_system_tag_attrib_seq_region_update_ids = []

        for seq_region_name_raw in chromosome_display_order_list:
            # wrap seq_region_name_raw into seq_region struct { "name": seq_region_name_raw } and
            # get seq_region_id (perhaps, by using unversioned name)
            seq_region_name, seq_region_id, unversioned_name = self.name_and_id_from_seq_region_item(
                {"name": seq_region_name_raw}, seq_region_map, try_unversion=unversion
            )

            # append trio to the resulting list
            regions_with_ranks_from_assembly_metadata.append(
                (seq_region_name, seq_region_id, unversioned_name)
            )

            # filling insert lists for "karyotype_rank" and "coord_system_tag" attributes
            rank_insertions_trios.append(
                (seq_region_id, karyotype_rank_attrib_id, len(rank_insertions_trios) + 1)
            )
            coord_system_tag_attrib_insertion_trios.append(
                (seq_region_id, coord_system_tag_attrib_id, self.quote_or_null(coord_system_tag))
            )

            # filling update list for "coord_system_tag" with seq_region_ids
            coord_system_tag_attrib_seq_region_update_ids.append(seq_region_id)

        # run insertion SQL for "karyotype_rank"
        self.insert_to_db(
            rank_insertions_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            self.pjc(work_dir, "karyotype_rank_insertion"),
            ignore=True,
        )

        # run insertion SQL for "coord_system_tag"
        self.insert_to_db(
            coord_system_tag_attrib_insertion_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            self.pjc(work_dir, "coord_system_tag_insertion"),
            ignore=True,
        )

        # forcing update of the "coord_system_tag"
        if force_update_coord_system_tag and coord_system_tag_attrib_seq_region_update_ids:
            seq_region_ids_str = ",".join(map(str, coord_system_tag_attrib_seq_region_update_ids))
            self.update_db_single_group(
                {"value": self.quote_or_null(coord_system_tag)},
                "seq_region_attrib",
                self.pjc(work_dir, "coord_system_tag_update"),
                where=f"attrib_type_id = {coord_system_tag_attrib_id} and seq_region_id in ({seq_region_ids_str})",
            )

        # return resulting list of regions with bands trios
        return regions_with_ranks_from_assembly_metadata

    def add_karyotype_rank_for_chromosomes(
        self, attrib_type_map: dict, work_dir: str, chromosome_coord_system_name="chromosome"
    ) -> list:  # [ (seq_region_name, seq_region_id, unversioned_name) ]
        """
        Add `karyotype_rank` attributes for seq region data from the "chromosome" coordinate system.

        Returns list of [ (seq_region_name, seq_region_id, unversioned_name) ] trios for seq_regions with updated karyotype_ranks.
        Not altering "coord_system_tag" tag attributes.

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """

        # resulting list of seq_region with karyotype_rank
        #   list of top level seq regions from the `chromosome_coord_system_name`
        chromomes_seq_regions = self.get_toplevel_from_cs(
            chromosome_coord_system_name, self.pjc(work_dir, "chromosome_seq_regions")
        )

        if not chromomes_seq_regions:
            return chromomes_seq_regions

        # technical / optimization. get external_db_id for "karyotype_rank"
        karyotype_rank_attrib_id = self.id_from_map_or_die(
            "karyotype_rank", attrib_type_map, "attrib_type_map"
        )

        # set/update proper attributes for "chromomosome" regions
        rank_insertions_trios = []
        for _, seq_region_id, _ in chromomes_seq_regions:
            rank_insertions_trios.append(
                (seq_region_id, karyotype_rank_attrib_id, len(rank_insertions_trios) + 1)
            )

        # run insertion SQL for "karyotype_rank"
        #    do not alter "coord_system_tag" anyhow in the case of the "chromosome" coord_system
        self.insert_to_db(
            rank_insertions_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            self.pjc(work_dir, "karyotype_rank_insertion"),
            ignore=True,
        )

        return chromomes_seq_regions

    def add_karyotype_rank_from_bands_info(
        self,
        regions_with_karyotype_bands: list,  # [ (seq_region_name, seq_region_id, unversioned_name) ]
        other_regions_with_ranks: list,  # [ (seq_region_name, seq_region_id, unversioned_name) ]
        attrib_type_map: dict,
        work_dir: str,
    ):
        """
        Add karyotype_ranks for `regions_with_karyotype_bands` (those with karyotype bands in seq_region metadata) but not present in `other_regions_with_ranks` list.

        Too close to the DB schema.
        """
        # form set of used seq_region_id(s)
        regions_with_ranks = frozenset(map(lambda el: el[1], other_regions_with_ranks))

        # get set of seq_region_ids with bands
        regions_with_bands = set(map(lambda el: el[1], regions_with_karyotype_bands))

        # seq_region_ids list to add ranks for
        region_ids_with_bands_but_no_karyotype_ranks = sorted(list(regions_with_bands - regions_with_ranks))

        # return if nothing to add
        if not region_ids_with_bands_but_no_karyotype_ranks:
            return

        # technical / optimization. get external_db_id for "karyotype_rank"
        karyotype_rank_attrib_id = self.id_from_map_or_die(
            "karyotype_rank", attrib_type_map, "attrib_type_map"
        )

        # set/update proper attributes for "chromomosome" regions
        rank_insertions_trios = []
        for seq_region_id in region_ids_with_bands_but_no_karyotype_ranks:
            rank_insertions_trios.append(
                (
                    seq_region_id,
                    karyotype_rank_attrib_id,
                    len(rank_insertions_trios) + 1 + len(regions_with_ranks),
                )
            )

        # run insertion SQL for "karyotype_rank"
        #    do not alter "coord_system_tag" anyhow in the case of the "chromosome" coord_system
        self.insert_to_db(
            rank_insertions_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            self.pjc(work_dir, "karyotype_rank_insertion"),
            ignore=True,
        )

    def unversion_scaffolds(self, cs_rank, logs):
        """
        Unversion scaffold, remove ".\d$" from seq_region.names if there's a need

        Non-versioned syns for contigs (lower, sequence level), versioned for the rest.
        """
        # coord_systems to ignore when adding synonyms for sequence_level cs
        artificial_cs = frozenset(self.param("artificial_cs").split(","))
        seq_cs, max_rank = max(
            [(c, r) for c, r in cs_rank.items() if c not in artificial_cs], key=lambda k: k[1]
        )
        for cs in cs_rank:
            if cs == seq_cs:
                # for non-sequence level cs, store original name (with version) as "sr_syn_src" synonym
                xdb = self.param("sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, self.pjc(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region_synonym", "synonym", self.pjc(logs, "unv_srs", cs))
            else:
                # for sequence-level cs, store original name (with version) as "versioned_sr_syn_src" synonym
                xdb = self.param("versioned_sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, self.pjc(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region", "name", self.pjc(logs, "unv_sr", cs))

    def coord_sys_order(self, cs_order_str):
        cs_order_lst = map(lambda x: x.strip(), cs_order_str.split(","))
        return {e: i for i, e in enumerate(filter(lambda x: len(x) > 0, cs_order_lst))}

    def used_cs_ranks(self, agps, cs_order, noagp_default=None):
        # rank cs_names, met in agps.keys ("-" separated), i.e. "scaffold-contig"
        #   only agps keys used, values are ignored
        #   use noagp_cs_name_default for "noagp" assemblies
        if agps is None:
            if noagp_default is None:
                raise Exception("NoAGP assembly with no default coordinate system name")
            cs_used_set = frozenset([noagp_default])
        else:
            cs_used_set = frozenset(sum(map(lambda x: x.split("-"), agps.keys()), []))

        cs_unknown = cs_used_set.difference(cs_order.keys())
        if len(cs_unknown) > 0:
            raise Exception("Unknown coordinate system(s) %s" % {str(cs_unknown)})
        return {e: i for i, e in enumerate(sorted(cs_used_set, key=lambda x: -cs_order[x]), start=1)}

    def chunk_contigs(self, fasta, cs_ranks, agps, work_dir, chunk_size=0, chunks_cs_name="ensembl_internal"):
        """
        chunk dna sequence fasta
          no chunking if chunk_size < 50k
        """
        chunk_size_min_len = self.param_required("sequence_data_chunk_min_len")
        if chunk_size < chunk_size_min_len:
            return fasta, cs_ranks, agps

        # split using script
        os.makedirs(work_dir, exist_ok=True)
        _stderr = f"{work_dir}/chunking.stderr"
        _out_agp = f"{work_dir}/chunks.agp"
        _out_fasta = f"{work_dir}/chunks.fasta"

        _splitter = r"fasta_chunk"
        split_cmd = f"{_splitter} --chunk_size {chunk_size} --agp_out {_out_agp} --out {_out_fasta} --fasta_dna {fasta} 2> {_stderr}"

        print(f"running {split_cmd}", file=sys.stderr)
        # NB throws CalledProcessError if failed
        sp.run(split_cmd, shell=True, check=True)

        # add rank for chunks
        _cs_name, _cs_rank = sorted(cs_ranks.items(), key=lambda k: k[1])[-1]
        cs_ranks[chunks_cs_name] = _cs_rank + 1

        # add agps entry
        if agps is None:
            agps = dict()
        agps[f"{_cs_name}-{chunks_cs_name}"] = _out_agp

        return _out_fasta, cs_ranks, agps

    def prune_agps(self, agps, cs_order, agps_pruned_dir, pruning=True):
        # when loading agp sort by:
        #   highest component (cmp) cs level (lowest rank)
        #   lowest difference between cs ranks (asm - cmp)
        #   i.e: chromosome-scaffold scaffold-chunk chromosome-chunk
        #   if no agps return empty pruned result

        if not agps:
            return None

        agp_levels_sorted = self.order_agp_levels(agps, cs_order)

        # prune agps
        agps_pruned = dict()
        used_components = set()
        if not pruning:
            used_components = None
        for asm_cmp in agp_levels_sorted:
            agp_file_src = agps[asm_cmp]
            agp_file_dst = self.pjc(agps_pruned_dir, asm_cmp + ".agp")
            if self.agp_prune(agp_file_src, agp_file_dst, used_components) > 0:
                agps_pruned[asm_cmp] = agp_file_dst
        return agps_pruned

    def order_agp_levels(self, agps, cs_order):
        # sort agp for loading by:
        #   highest component (cmp) cs level (lowest rank)
        #   lowest difference between cs ranks (asm - cmp)
        #   i.e: chromosome-scaffold scaffold-chunk chromosome-chunk
        if not agps:
            return []

        agp_cs_pairs = list(map(lambda x: [x] + x.split("-"), agps.keys()))
        agp_levels = [(x[0], cs_order[x[1]], cs_order[x[2]]) for x in agp_cs_pairs]

        bad_agps = list(filter(lambda x: x[1] < x[2], agp_levels))
        if len(bad_agps) > 0:
            raise Exception("component cs has higher order than assembled cs %s" % (str(bad_agps)))

        agp_levels_sorted = [e[0] for e in sorted(agp_levels, key=lambda x: (-x[2], x[1] - x[2]))]
        return agp_levels_sorted

    def load_seq_data(self, fasta, agps, cs_rank, log_pfx):
        """loads sequence data for various coordinate systems accordingly with their rank"""
        asm_v = self.asm_name()

        sequence_rank = max(cs_rank.values())
        for cs, rank in sorted(cs_rank.items(), key=lambda p: -p[1]):
            logs = self.pjc(log_pfx, "%02d_%s" % (rank, cs))
            if rank == sequence_rank:
                self.load_cs_data(cs, rank, "fasta", asm_v, fasta, logs, loaded_regions=None, seq_level=True)
            else:
                useful_agps = list(filter(lambda x: cs in x, agps and agps.keys() or []))
                if len(useful_agps) == 0:
                    raise Exception("non-seq_level cs %s has no agps to assemble it from" % (cs))
                loaded_regions = set()
                for pair, agp_file_pruned in map(lambda k: (k, agps[k]), useful_agps):
                    if not pair.startswith(cs + "-"):
                        continue
                    self.load_cs_data(cs, rank, pair, asm_v, agp_file_pruned, logs, loaded_regions)

    def load_cs_data(self, cs, rank, pair, asm_v, src_file, log_pfx, loaded_regions=None, seq_level=False):
        """creates a coord_system and loads sequence or assembly(AGP) data for corresponding seqregions

        doesn't load already seen sequences
        """
        # NB load_seq_region.pl and load_agp.pl are not failing on parameter errors (0 exit code)
        os.makedirs(dirname(log_pfx), exist_ok=True)
        additional_load = self.param_bool("load_additional_sequences")
        if seq_level:
            self.load_seq_region(cs, rank, asm_v, src_file, log_pfx, seq_level, additional_load)
        elif loaded_regions is not None:
            new_regions = set()
            clean_file = src_file + ".regions_deduped"
            self.filter_already_loaded_regions_from_agp(src_file, clean_file, loaded_regions, new_regions)
            self.load_seq_region(cs, rank, asm_v, clean_file, log_pfx, seq_level, additional_load)
            loaded_regions.update(new_regions)
        if not seq_level:
            self.load_agp(pair, asm_v, src_file, log_pfx)

    def filter_already_loaded_regions_from_agp(self, src_file, dst_file, loaded_regions, new_regions):
        with open(src_file) as src:
            with open(dst_file, "w") as dst:
                for line in src:
                    fields = line.strip().split("\t")
                    (
                        asm_id,
                        asm_start,
                        asm_end,
                        asm_part,
                        type_,
                        cmp_id,
                        cmp_start,
                        cmp_end,
                        cmp_strand,
                    ) = fields
                    if type_ in "NU" or asm_id in loaded_regions:
                        continue
                    new_regions.add(asm_id)
                    print(line.strip(), file=dst)

    def agp_prune(self, from_file: str, to_file: str, used: set = None):
        """
        Remove already components from the AGP file if they are seen in "used" set
        """
        # reomve used component
        #   and GAPS as they are not used by 'ensembl-analysis/scripts/assembly_loading/load_agp.pl'
        os.makedirs(dirname(to_file), exist_ok=True)
        open_ = self.is_gz(from_file) and gzip.open or open
        if used is None:
            cmd = r"""{_cat} {_file} > {_out}""".format(
                _cat=self.is_gz(from_file) and "zcat" or "cat", _file=from_file, _out=to_file
            )
            print("running %s" % (cmd), file=sys.stderr)
            sp.run(cmd, shell=True, check=True)
            return 1
        writes = 0
        with open_(from_file, "r") as src:
            with open(to_file, "w") as dst:
                for line in src:
                    fields = line.strip().split("\t")
                    (
                        asm_id,
                        asm_start,
                        asm_end,
                        asm_part,
                        type_,
                        cmp_id,
                        cmp_start,
                        cmp_end,
                        cmp_strand,
                    ) = fields
                    if type_ in "NU" or cmp_id in used:
                        continue
                    used.add(cmp_id)
                    print(line.strip(), file=dst)
                    writes += 1
        return writes

    def get_external_db_mapping(self) -> dict:
        """
        Get a map from a file for external_dbs to Ensembl dbnames from "external_db_map" module(!) param
        """
        external_map_path = self.param("external_db_map")
        db_map = dict()
        if external_map_path is None:
            return db_map

        # Load the map
        with open(external_map_path, "r") as map_file:
            for line in map_file:
                if line.startswith("#"):
                    continue
                line = re.sub(r"#.*", "", line)
                if re.match(r"^\s*$", line):
                    continue
                (from_name, to_name, *rest) = line.strip().split("\t")
                if len(rest) > 0 and rest[0].upper() != "SEQ_REGION":
                    continue
                if to_name == "_IGNORE_":
                    continue
                db_map[from_name] = to_name
        return db_map

    # UTILS
    def db_string(self):
        return "-dbhost {host_} -dbport {port_} -dbuser {user_} -dbpass {pass_} -dbname {dbname_} ".format(
            host_=self.param("dbsrv_host"),
            port_=self.param("dbsrv_port"),
            user_=self.param("dbsrv_user"),
            pass_=self.param("dbsrv_pass"),
            dbname_=self.param("db_name"),
        )

    def pjc(self, *parts: list) -> str:
        """
        Join path parts and try to create every directory but the last one.
        """
        if not parts:
            return None

        parts = list(parts)
        last = parts.pop()
        prefix = pj(*parts)

        os.makedirs(prefix, exist_ok=True)

        return pj(prefix, last)

    def is_gz(self, filename):
        return filename.endswith(".gz")

    def asm_name(self):
        asm = self.from_param("genome_data", "assembly")
        if "name" not in asm:
            raise Exception("no assembly/name in genome_data")
        return asm["name"]

    # TODO: add some metafunc setter getter
    def from_param(self, param, key, not_throw=False):
        data = self.param_required(param)
        if key not in data:
            if not_throw:
                return None
            else:
                raise Exception("Missing required %s data: %s" % (param, key))
        return data[key]

    def param_bool(self, param):
        val = self.param(param)
        return bool(val) and "0" != val

    def load_map_from_sql_stdout(self, in_file, skip_header=False):
        """
        Load map from the SQL output

        Process input in_file with "key  value" pairs and load then
        into the {key : value} map.
        Skips header if skip_header.
        """
        data = dict()
        with open(in_file) as pairs_file:
            for line in pairs_file:
                if skip_header:
                    skip_header = False
                    continue
                (key, val) = line.strip().split("\t")
                data[key] = val
        return data

    def name_and_id_from_seq_region_item(
        self,
        seq_region_item: dict,
        seq_region_map: dict,
        try_unversion: bool = False,
        throw_missing: bool = True,
    ) -> (str, str, str):
        """
        Get (seq_region_name, seq_region_id, unversioned_name) from seq_region_item struct(dict)

        Gets unversioned_name only if "try_unversion" is True.
        Throws exception if not able to get seq_region_id from "seq_region_map" and "throw_missing" is true.
        """
        #   get seq_region_id (perhaps, by using unversioned name)
        seq_region_name = seq_region_item["name"]
        seq_region_id = seq_region_map.get(seq_region_name, None)
        unversioned_name = None
        if seq_region_id is None and try_unversion:
            # try to get seq_region_id for the unversioned name
            unversioned_name = re.sub(r"\.\d+$", "", seq_region_name)
            seq_region_id = seq_region_map.get(unversioned_name, "")

        # oops, we don't know such seq_region name
        if not seq_region_id and throw_missing:
            raise Exception(f"Not able to find seq_region for '{seq_region_name}'")

        return (seq_region_name, seq_region_id, unversioned_name)

    def id_from_map_or_die(self, key: str, map_dict: dict, name_for_panic):
        value = map_dict.get(key, None)
        if value is None:
            raise Exception(f"no such key '{key}' in '{name_for_panic}' map")
        return value

    ## Utilities using external scripts
    def remove_IUPAC(self, from_file: str, to_file: str):
        """remove non-valid symbols from FASTA file (using sed) ans store the result in a different location"""
        IUPAC = self.param("IUPAC")
        os.makedirs(dirname(to_file), exist_ok=True)
        cmd = r"""{_cat} {_file} | sed -r '/^[^>]/ {{ s/[{_IUPAC}]/N/g; s/{_iupac}/n/g }}' > {_out}""".format(
            _cat=self.is_gz(from_file) and "zcat" or "cat",
            _file=from_file,
            _IUPAC=IUPAC.upper(),
            _iupac=IUPAC.lower(),
            _out=to_file,
        )
        print("running %s" % (cmd), file=sys.stderr)
        return sp.run(cmd, shell=True, check=True)

    def load_seq_region(
        self,
        cs: str,
        rank: str,
        asm_v: str,
        src_file: str,
        log_pfx: str,
        seq_level=False,
        additional_load=False,
    ):
        """ensembl-analysis script (load_seq_region.pl) based utility for loading seq_regions FASTA sequences"""
        en_root = self.param_required("ensembl_root_dir")
        cmd = (
            r"""{_loader} {_db_string} {_asm_v_flag} -default_version -ignore_ambiguous_bases """
            + r"""    -rank {_rank} -coord_system_name {_cs} {_sl_flag} -{_tag}_file {_file}"""
            + r"""     > {_log}.stdout 2> {_log}.stderr"""
        ).format(
            _loader="perl %s"
            % (self.pjc(en_root, r"ensembl-analysis/scripts/assembly_loading/load_seq_region.pl")),
            _db_string=self.db_string(),
            _asm_v_flag=not additional_load and f"-coord_system_version {asm_v}" or "",
            _rank=rank,
            _cs=cs,
            _sl_flag=seq_level and "-sequence_level" or "",
            _tag=seq_level and "fasta" or "agp",
            _file=src_file,
            _log="%s_seq" % (log_pfx),
        )
        print("running %s" % (cmd), file=sys.stderr)
        return sp.run(cmd, shell=True, check=True)

    def load_agp(self, pair, asm_v, src_file, log_pfx):
        """ensembl script (load_agp.pl) based utility for loading seq_regions assembly data (AGPs)"""
        en_root = self.param_required("ensembl_root_dir")
        (asm_n, cmp_n) = pair.strip().split("-")
        cmd = (
            r"""{_loader} {_db_string} -assembled_version {_asm_v} """
            + r"""    -assembled_name {_asm} -component_name {_cmp} """
            + r"""    -agp_file {_file} """
            + r"""    > {_log}.stdout 2> {_log}.stderr"""
        ).format(
            _loader="perl %s" % (self.pjc(en_root, r"ensembl-analysis/scripts/assembly_loading/load_agp.pl")),
            _db_string=self.db_string(),
            _asm_v=asm_v,
            _asm=asm_n,
            _cmp=cmp_n,
            _file=src_file,
            _log="%s_agp_%s" % (log_pfx, pair.replace("-", "_")),
        )
        print("running %s" % (cmd), file=sys.stderr)
        return sp.run(cmd, shell=True, check=True)

    def set_toplevel(self, log_pfx, ignored_cs=[]):
        """
        Set toplevel(6) seq_region_attrib using ensembl script.

        Uses set_toplevel.pl ensembl script.
        """
        # set top_level(6) seq_region_attrib
        os.makedirs(dirname(log_pfx), exist_ok=True)
        en_root = self.param_required("ensembl_root_dir")
        cmd = (
            r"""{_set_tl} {_db_string} {_ignored_cs} """ + r"""     > {_log}.stdout 2> {_log}.stderr"""
        ).format(
            _set_tl="perl %s"
            % (self.pjc(en_root, r"ensembl-analysis/scripts/assembly_loading/set_toplevel.pl")),
            _db_string=self.db_string(),
            _ignored_cs=" ".join(map(lambda x: "-ignore_coord_system %s" % (x), ignored_cs)),
            _log=log_pfx,
        )
        print("running %s" % (cmd), file=sys.stderr)
        sp.run(cmd, shell=True, check=True)

        # remove toplevel attribute for seq_regions that are components
        self.remove_components_from_toplevel(log_pfx)

    ## SQL executor and utilities using plain SQL
    def run_sql_req(self, sql, log_pfx, from_file=False):
        os.makedirs(dirname(log_pfx), exist_ok=True)

        sql_option = r""" -sql '{_sql}' """.format(_sql=sql)
        if from_file:
            sql_option = r""" < '{_sql}' """.format(_sql=sql)

        cmd = r"""{_dbcmd} -url "{_srv}{_dbname}" {_sql_option} > {_out} 2> {_err}""".format(
            _dbcmd="perl %s/scripts/db_cmd.pl" % os.getenv("EHIVE_ROOT_DIR"),
            _srv=self.param("dbsrv_url"),
            _dbname=self.param("db_name"),
            _sql_option=sql_option,
            _out=log_pfx + ".stdout",
            _err=log_pfx + ".stderr",
        )
        print("running %s" % (cmd), file=sys.stderr)
        return sp.run(cmd, shell=True, check=True)

    def add_contig_ena_attrib(self, log_pfx, cs_name="contig"):
        """
        Add ENA attrib for contigs if their names are ENA accessions

        Nno sequence_level checks are used -- just cs name.
        See ensembl-datacheck/lib/Bio/EnsEMBL/DataCheck/Checks/SeqRegionNamesINSDC.pm .
        SQL code.
        """
        sql = r"""insert ignore into seq_region_attrib (seq_region_id, attrib_type_id, value)
                select
                  sr.seq_region_id, at.attrib_type_id, "ENA"
                from
                  seq_region sr, coord_system cs, attrib_type at
                where   sr.coord_system_id = cs.coord_system_id
                    and cs.name = "%s"
                    and at.code = "external_db"
              ;""" % (
            cs_name
        )
        return self.run_sql_req(sql, log_pfx)

    def copy_sr_name_to_syn(self, cs, x_db, log_pfx):
        """
        Store original seq_region names as seq_region_synonym

        Store original seq_region names (from a given cood_systen, "cs" param) as seq_region_synonyms (using "x_db" external source name)
        SQL code.
        """
        asm_v = self.asm_name()
        sql = r"""insert into seq_region_synonym (seq_region_id, synonym, external_db_id)
                  select
                      sr.seq_region_id, sr.name, xdb.external_db_id
                  from
                     seq_region sr, external_db xdb, coord_system cs
                  where   xdb.db_name = "%s"
                      and sr.coord_system_id = cs.coord_system_id
                      and cs.name = "%s"
                      and cs.version = "%s"
                      and sr.name like "%%._"
                ;""" % (
            x_db,
            cs,
            asm_v,
        )
        return self.run_sql_req(sql, log_pfx)

    def add_asm_mappings(self, cs_pairs, log_pfx):
        """
        Adds "assembly.mapping" strings to meta table.

        Nullifies asm_mappings contig versions as well, but don't nullify toplevel
        Doesn't add mapping id there is a single CS
        SQL code.
        """
        if cs_pairs is None or len(cs_pairs) < 1:
            return
        asm_v = self.asm_name()
        for pair in cs_pairs:
            higher, lower = pair.strip().split("-")
            sql = r"""insert ignore into meta (species_id, meta_key, meta_value) values
                    (1, "assembly.mapping", "{_higher}:{_v}|{_lower}:{_v}")
                  ;""".format(
                _v=asm_v, _higher=higher, _lower=lower
            )
            self.run_sql_req(sql, self.pjc(log_pfx, pair))

    def remove_components_from_toplevel(self, log_pfx):
        """
        Remove toplevel attribute for seq_regions that are "components" (parts of different seq_regions).

        SQL code.
        """

        # get list of seq_regions that are components
        sql_not_toplevel_list = r"""
          select distinct a.cmp_seq_region_id, sr_c.name, cs_c.name
            from assembly a,
                 seq_region sr_c, seq_region sr_a,
                 coord_system cs_c, coord_system cs_a
            where a.cmp_seq_region_id = sr_c.seq_region_id
              and a.asm_seq_region_id = sr_a.seq_region_id
              and sr_c.coord_system_id = cs_c.coord_system_id
              and sr_a.coord_system_id = cs_a.coord_system_id
              and cs_c.attrib like "%default_version%"
              and cs_a.attrib like "%default_version%"
              and sr_c.seq_region_id in (
                select distinct seq_region_id
                  from seq_region_attrib
                  where attrib_type_id = 6 and value = 1
              )
        """
        # perhaps, make sense to check sr_c.coord_system_id != sr_a.coord_system_id
        self.run_sql_req(sql_not_toplevel_list, ".".join([log_pfx, "not_toplevel_list"]))

        # delete wrongly assigned attribs
        sql_not_toplevel_delete = r"""
          delete from seq_region_attrib
            where attrib_type_id = 6 and value = 1
              and seq_region_id in (
                select distinct a.cmp_seq_region_id
                  from assembly a,
                       seq_region sr_c, seq_region sr_a,
                       coord_system cs_c, coord_system cs_a
                  where a.cmp_seq_region_id = sr_c.seq_region_id
                    and a.asm_seq_region_id = sr_a.seq_region_id
                    and sr_c.coord_system_id = cs_c.coord_system_id
                    and sr_a.coord_system_id = cs_a.coord_system_id
                    and cs_c.attrib like "%default_version%"
                    and cs_a.attrib like "%default_version%"
              );
        """
        # perhaps, check sr_c.coord_system_id != sr_a.coord_system_id
        self.run_sql_req(sql_not_toplevel_delete, ".".join([log_pfx, "not_toplevel_delete"]))

    def sr_name_unversion(self, cs, tbl, fld, log_pfx):
        """
        Remove version suffix from the seq_region names

        Removes '\.\d+$' suffices from the seq_region names
        SQL code.
        """
        # select synonym, substr(synonym,  1, locate(".", synonym, length(synonym)-2)-1)
        #     from seq_region_synonym  where synonym like "%._"
        asm_v = self.asm_name()
        sql = r"""update {_tbl} t, seq_region sr, coord_system cs
                    set
                      t.{_fld} = substr(t.{_fld},  1, locate(".", t.{_fld}, length(t.{_fld})-2)-1)
                    where t.{_fld} like "%._"
                      and t.seq_region_id = sr.seq_region_id
                      and sr.coord_system_id = cs.coord_system_id
                      and cs.name = "{_cs}"
                      and cs.version = "{_asm_v}"
                ;""".format(
            _tbl=tbl, _fld=fld, _cs=cs, _asm_v=asm_v
        )
        return self.run_sql_req(sql, log_pfx)

    def nullify_ctg_cs_version(self, cs_order, log_pfx: str):
        """
        Nullify every CS version with rank larger than that of "contig", but don't nullify toplevel ones.

        SQL code
        """
        asm_v = self.asm_name()
        # get cs_info (and if they have toplevel regions)
        sql = r"""select cs.coord_system_id as coord_system_id,
                         cs.name, cs.rank, (tl.coord_system_id is NULL) as no_toplevel
                    from coord_system cs
                      left join (
                        select distinct sr.coord_system_id
                          from seq_region sr, seq_region_attrib sra, attrib_type at
                          where at.code = "toplevel"
                            and sra.attrib_type_id = at.attrib_type_id
                            and sra.value = 1
                            and sra.seq_region_id = sr.seq_region_id
                      ) as tl on tl.coord_system_id = cs.coord_system_id
                    where cs.version = "{_asm_v}"
                    order by rank
              ;""".format(
            _asm_v=asm_v
        )
        # run_sql
        toplvl_pfx = self.pjc(log_pfx, "toplvl_info")
        self.run_sql_req(sql, toplvl_pfx)
        # load info
        cs_info = []
        with open(toplvl_pfx + ".stdout") as f:
            header = None
            for line in f:
                if header is None:
                    header = line.strip().split("\t")
                    continue
                cs_info.append(dict(zip(header, line.strip().split())))
        # return if there's no coord_systems to nullify versions for
        if not cs_info:
            return
        # get list of known cs from cs_order to clean version from
        nullify_cs_version_from = self.param("nullify_cs_version_from")
        if not nullify_cs_version_from or nullify_cs_version_from not in cs_order:
            return
        cs_thr_index = cs_order[nullify_cs_version_from]
        cs_names_to_keep_ver = frozenset([nm for (nm, ind) in cs_order.items() if ind > cs_thr_index])

        # choose cs rank threshold to start clearing version from
        clear_lst = [
            (cs["coord_system_id"], cs["name"])
            for cs in cs_info
            if (bool(int(cs["no_toplevel"])) and cs["name"] not in cs_names_to_keep_ver)
        ]

        # run sql
        if clear_lst:
            clear_pfx = self.pjc(log_pfx, "clear")
            with open(clear_pfx + ".sql", "w") as clear_sql:
                for cs_id, cs_name in clear_lst:
                    sql = r"""
                        update meta set
                            meta_value=replace(meta_value, "|{_cs_name}:{_asm_v}", "|{_cs_name}")
                            where meta_key="assembly.mapping";
                        update coord_system set version = NULL where coord_system_id = {_cs_id};
                    """.format(
                        _asm_v=asm_v, _cs_name=cs_name, _cs_id=cs_id
                    )
                    print(sql, file=clear_sql)
            self.run_sql_req(clear_pfx + ".sql", clear_pfx, from_file=True)

    def load_map_from_core_db(self, table, cols, work_dir) -> dict:
        """
        Load 2 "cols" from core db "table" as map

        Load { cols[0] : cols[1] } map from the core db "table"
        SQL code
        """
        out_pfx = self.pjc(work_dir, f"{table}_map")
        sql = f"""select {cols[0]}, {cols[1]} FROM {table};"""
        res = self.run_sql_req(sql, out_pfx)

        out_file = out_pfx + ".stdout"
        data = self.load_map_from_sql_stdout(out_file, skip_header=True)
        if not data:
            raise Exception(f"No '{table}' map loaded from '{out_file}'")
        return data

    def load_seq_region_synonyms_trios_from_core_db(self, work_dir: str) -> list:
        # was get_db_syns
        """
        Load seq_region_synonyms from from core db into [(seq_region_id, name, synonym)...] list

        SQL code
        """
        out_pfx = self.pjc(work_dir, f"seq_region_synonyms")
        sql = r"""select sr.seq_region_id as seq_region_id, sr.name, srs.synonym
                 from seq_region sr left join seq_region_synonym srs
                 on sr.seq_region_id = srs.seq_region_id
                 order by sr.seq_region_id
              ;"""

        res = self.run_sql_req(sql, out_pfx)

        syn_trios = []
        out_file = out_pfx + ".stdout"
        with open(out_file) as syns_file:
            skip_header = True
            for line in syns_file:
                if skip_header:
                    skip_header = False
                    continue
                (sr_id, name, syn) = line.strip().split("\t")
                syn_trios.append((sr_id, name, syn))
        return syn_trios

    def insert_to_db(
        self, list_of_tuples: list, table_name: str, col_names: list, work_dir: str, ignore: bool = True
    ):
        """
        Insert into the core db's {table_name} tuples from {list_of_tuples} as col_names.

        Use `quote_or_null` (see definition below) method for string values, when putting values into `list_of_tuples`
        SQL code
        """
        # return if nothing to do
        if not list_of_tuples:
            return

        # prepare request parts
        ignore_str = ignore and "IGNORE" or ""
        cols_str = ", ".join(col_names)

        # generate file with the insert SQL command
        insert_sql_file = self.pjc(work_dir, "insert.sql")
        with open(insert_sql_file, "w") as sql:
            print(f"INSERT {ignore_str} INTO {table_name} ({cols_str}) VALUES", file=sql)
            values_sep = ""
            for tpl in list_of_tuples:
                tpl_str = ", ".join(map(str, tpl))
                print(f"{values_sep}({tpl_str})", file=sql)
                values_sep = ", "
            print(";", file=sql)

        # run insert SQL from file
        self.run_sql_req(insert_sql_file, self.pjc(work_dir, "insert"), from_file=True)

    def quote_or_null(self, val: str, quotes: str = "'", null: str = "NULL", strings_only=True) -> str:
        """
        Return `val` wrapped in `quotes` or `null` value

        Quotes only strings (instances of `str`) if strings_only is True.
        """
        if val is None:
            return null
        if strings_only and isinstance(val, str):
            return f"{quotes}{val}{quotes}"
        return val

    def update_db_single_group(
        self, dict_of_col_to_value: dict, table_name: str, work_dir: str, where: str = None
    ):
        """
        Update given `table` name in db; set `col = val` for all key/value pairs from `dict_of_cols_to_values`

        If `where` condition is present its value is used for the "WHERE" SQL clause.
        Use `quote_or_null` (see definition below) method for string values, when putting values into `list_of_tuples`

        SQL code
        """
        # return if nothing to do
        if not dict_of_col_to_value:
            return

        # prepare request parts
        where_str = where and f"WHERE {where}" or ""
        col_val_str = ", ".join([f"{col} = {val}" for col, val in dict_of_col_to_value.items()])

        # generate file with the insert SQL command
        update_sql_file = self.pjc(work_dir, "update.sql")
        with open(insert_sql_file, "w") as sql:
            print(f"UPDATE {table_name} SET {col_val_str} {where_str};", file=sql)
            values_sep = ""
            for tpl in list_of_tuples:
                tpl_str = ", ".join(map(str, tpl))
                print(f"{values_sep}({tpl_str})", file=sql)
                values_sep = ", "
            print(";", file=sql)

        # run insert SQL from file
        self.run_sql_req(update_sql_file, self.pjc(work_dir, "update"), from_file=True)

    def get_toplevel_from_cs(self, coord_system_name, work_dir) -> list:
        """
        Returns list of [ (seq_region_name, seq_region_id, "") ] trios for toplevel seq_regions from coord system with `coord_system_name`
          or having  "coord_system_tag" attribute with the `coord_system_name` value

        SQL code
        """
        out_pfx = self.pjc(work_dir, f"toplevel_from_{coord_system_name}")
        sql = f"""SELECT DISTINCT sr.name, sr.seq_region_id
                FROM seq_region sr,
                     seq_region_attrib sra,
                     coord_system cs,
                     attrib_type at
                WHERE sr.seq_region_id = sra.seq_region_id
                  AND sr.coord_system_id = cs.coord_system_id
                  AND sra.attrib_type_id = at.attrib_type_id
                  AND (  ( cs.name = "{coord_system_name}" and at.code = "toplevel" )
                      OR ( at.code = "coord_system_tag" and sra.value = "{coord_system_name}" )
                      )
                  ORDER BY sr.seq_region_id;
               """

        res = self.run_sql_req(sql, out_pfx)

        sr_trios = []
        out_file = out_pfx + ".stdout"
        with open(out_file) as sr_file:
            skip_header = True
            for line in sr_file:
                if skip_header:
                    skip_header = False
                    continue
                (name, sr_id) = line.strip().split("\t")
                sr_trios.append((name, sr_id, ""))

        return sr_trios
