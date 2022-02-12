#!env python3

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
            'cs_order' : 'chunk,contig,supercontig,non_ref_scaffold,scaffold,primary_assembly,superscaffold,linkage_group,chromosome',

            'IUPAC' : 'RYKMSWBDHV',    # symbols to be replaced with N in the DNA sequences (ensembl core(107) doesn't support the whole IUPAC alphabet for DNA)

            # unversion scaffold, remove ".\d$" from seq_region.names if there's a need
            'unversion_scaffolds' : 0,
            'versioned_sr_syn_src' : 'INSDC', # INSDC(50710) # if unversioning non-sequence level cs, store original name (with version) as this synonym
            'sr_syn_src' : 'BRC4_Community_Symbol', # BRC4_Community_Symbol(211) # if unversioning sequence-level cs, store original name (with version) as this synonym

            # nullify coord_system version for the given coord_system name 
            'nullify_cs_version_from' : 'contig',

            # default coord_system name for single-level (no AGPs assemblies)
            'noagp_cs_name_default' : 'primary_assembly',

            # file for additional mapping of the synonym sources to the external_db (ensembl), as used by "get_external_db_mapping" function below
            'external_db_map' : None,

            # set "coord_system_tag" seq_region attribute to this value if there's a corresponding "chromosome_display_order" list in genome.json metadata
            #   if None, only "chromosome" coord system is processed (if present)
            #   (see add_chr_karyotype_rank definition below )
            'cs_tag_for_ordered' : None,

            # BRC4 compatibility mode; if on, "(EBI|BRC4)_seq_region_name" seq_region_attributes are added.
            #   Blocked by the "swap_gcf_gca" option. In this case insertion should be done on later pipeline stage after seq_region name swapping.
            'brc4_mode' : True,

            # Whether to use RefSeq names as additional seq_region synonyms (if available) or not (see add_sr_synonyms definition below)
            #  Does not swap anything actually, just loads synonyms to be used by a later "swapping" stage
            #  Disables BRC4 compatibilty mode (see the "brc4_mode" option comment).
            'swap_gcf_gca' : False,

            # list of coord systems used in "-ignore_coord_system" options of the "ensembl-analysis/scripts/assembly_loading/set_toplevel.pl" script 
            #   part of the loading process
            'not_toplevel_cs' : [], # i.e. "contig", "non_ref_scaffold"

            # explicit list of seq_region properties (keys) to load as seq_region_attrib s (values) (see add_sr_attribs definition below) 
            #   if a dict's used as a value, treat its keys as "json_path" (/ as delim) map, i.e.
            #       { "added_sequence" : { "assembly_provider" : { "name" : ... } } } -> "added_sequence/assembly_provider/name"
            #   only flattable properties can be used, no arrays
            #   arrays should be processed separately (see `add_sr_synonyms` or `add_karyotype_bands` definitions)
            # see schema/seq_region_schema.json
            'sr_attrib_types' : {
                'circular' : 'circular_seq',
                'codon_table' : 'codon_table',
                'location' : 'sequence_location',
                'non_ref' : 'non_ref',
                'coord_system_level' : 'coord_system_tag',
                'added_sequence' : {
                    # json_path to attrib_type_code(str) mapping
                    'added_sequence/accession' :  'added_seq_accession',
                    'added_sequence/assembly_provider/name' :  'added_seq_asm_pr_nam',
                    'added_sequence/assembly_provider/url' :  'added_seq_asm_pr_url',
                    'added_sequence/annotation_provider/name' :  'added_seq_ann_pr_nam',
                    'added_sequence/annotation_provider/url' :  'added_seq_ann_pr_url',
                }, # added_sequence
            },
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
        external_db_map = self.load_map_from_core_db("external_db", ["db_name", "external_db_id"], work_dir) # for external_db
        attrib_type_map = self.load_map_from_core_db("attrib_type",  ["code", "attrib_type_id"], work_dir) # for attrib_type
        seq_region_map = self.load_map_from_core_db("seq_region",  ["name", "seq_region_id"], work_dir) # for seq_region

        # update synonyms and seq_region_attribs
        unversion = self.param("unversion_scaffolds")
        is_primary_assembly = self.from_param("manifest_data", "agp", not_throw = True) is None
        seq_region_file = self.from_param("manifest_data", "seq_region", not_throw = True)

        #   add seq_region synonyms
        self.add_sr_synonyms(seq_region_file,
                             seq_region_map,
                             external_db_map,
                             pj(work_dir, "seq_region_syns"),
                             unversion = unversion)

        #   add seq_region attributes
        self.add_sr_attribs(seq_region_file,
                            seq_region_map,
                            attrib_type_map,
                            pj(work_dir, "seq_region_attr"),
                            unversion = unversion)

        #   add seq_region EBI and BRC4 name attributes in the "BRC4 mode"
        #     special case of attributes adding with default values derived from seq_region names
        #     do not add if preparing to swap RefSeq and GeneBank ids; in this case attributes to be added at a later stage in pipeline
        #     (easier to insert then to update)
        if self.param("brc4_mode") and not self.param("swap_gcf_gca"):
            self.add_sr_ebi_brc4_names(seq_region_file,
                                       seq_region_map,
                                       attrib_type_map,
                                       pj(work_dir, "seq_region_ebi_brc4_name"),
                                       unversion = unversion)

        # add karyotype related data
        self.add_karyotype_data(seq_region_file,
                                seq_region_map,
                                attrib_type_map,
                                pj(work_dir, "karyotype"),
                                unversion = unversion)


    def initial_sequence_loading(self, work_dir: str):
        """
        initial preparation and loading of AGPs and fasta data

        initial preparation and loading of AGPs and fasta data using ensembl-analysis perl scripts
        """
        # preprocess FASTA with sequences
        #   rename IUPAC to N symbols using sed
        fasta_raw = self.from_param("manifest_data", "fasta_dna")
        fasta_clean = pj(work_dir, "fasta", "seq_no_iupac.fasta")
        self.remove_IUPAC(fasta_raw, fasta_clean)

        # start coord system ranking and agps processing
        agps = self.from_param("manifest_data", "agp", not_throw = True)

        # rank cs_names, met in agps.keys ("-" separated, i.e. "scaffold-contig") based on cs_order
        #   use noagp_cs_name_default for "noagp" assemblies
        cs_order = self.coord_sys_order(self.param("cs_order"))
        noagps_cs = self.param("noagp_cs_name_default")
        cs_rank = self.used_cs_ranks(agps, cs_order, noagps_cs)

        # remove gaps and lower_level mappings if the are coveres by higher level ones
        #   i.e.: remove 'contigN to chromosomeZ', if 'contigN to scaffoldM' and 'scaffoldM to chromosomeZ' are in place
        #   returns None if no agps provided
        agps_pruned_dir = pj(work_dir, "agps_pruned")
        agps_pruned = self.prune_agps(agps, cs_order, agps_pruned_dir, self.param_bool("prune_agp"))

        # empty agps_pruned ignored
        self.load_seq_data(fasta_clean, agps_pruned, cs_rank, pj(work_dir, "load"))

        # mark all the "contig"s or noagp_cs as being sourced from ENA
        if not self.param_bool("no_contig_ena_attrib"):
            if agps is None:
                self.add_contig_ena_attrib(pj(work_dir, "load", "set_ena"), cs_name = noagps_cs)
            else:
                self.add_contig_ena_attrib(pj(work_dir, "load", "set_ena"))

        # unversion scaffold, remove ".\d$" from names if there's a need
        if self.param_bool("unversion_scaffolds"):
            self.unversion_scaffolds(cs_rank, pj(work_dir, "unversion_scaffolds"))

        # add assembly mappings between various cs to meta table for the mapper to work properly
        cs_pairs = agps_pruned and agps_pruned.keys() or None
        self.add_asm_mappings(cs_pairs, pj(work_dir, "asm_mappings"))

        # set toplevel seq_region attribute
        self.set_toplevel(pj(work_dir, "set_toplevel"), self.param("not_toplevel_cs"))

        # nullify contig version and update mappings strings accordingly
        self.nullify_ctg_cs_version(pj(work_dir, "asm_mapping", "nullify_cs_versions"))


    def add_sr_synonyms(self,
                        seq_region_file: str,
                        seq_region_map: dict,
                        external_db_map: dict,
                        work_dir: str,
                        unversion: bool = False,
                        unversionable_sources_set: frozenset = frozenset(["INSDC", "RefSeq"])):
        """
        Add seq_region_synonym from the seq_region_file meta data file.

        Add seq_region_synonym from the schema/seq_region_schema.json compatible meta data file.
        Merge with the already exinsting ones in the db.
        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible
          * the unversioned synonyms from the unversionable_sources_set will be added as well as the original ones

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file: return

        # get seq_region ids, names, syns from db
        synonyms_trios_db = self.load_seq_region_synonyms_trios_from_core_db(pj(work_dir, "syns_from_core"))

        # form set of synonyms already present in db
        synonyms_in_db = frozenset( [trio[2] for trio in synonyms_trios_db if trio[2] != "NULL"] )

        # subset of sources to use the allowed the unversion synonyms
        unversionable_sources = unversion and unversionable_sources_set or frozenset()

        # technical / optimization. get external_db_id for "ensembl_internal_synonym"
        ensembl_internal_synonym_ext_db_id = self.id_from_map_or_die("ensembl_internal_synonym", external_db_map, "external_db_map")

        # get dict for additional mapping for sources to external_db names if there's one specified by "external_db_map" module param
        #   not to be confused with the external_db_map function parameter above
        additional_sources_mapping = self.get_external_db_mapping()

        # load synonyms from the json file
        synonyms_from_json = [] # [ (seq_region_id, synonym, external_db_id)... ] list of trios for inserting into db 
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in filter(lambda sr: "synonyms" in sr, seq_regions):
                # iterate through all seq_regions having "synonyms" 
                seq_region_name, seq_region_id, _ = \
                    self.name_and_id_from_seq_region_item(seq_region, try_unversion = unversion)

                # fill synonyms_from_json list of trios
                for synonym_item in list(seq_region["synonyms"]):
                    synonym_name = synonym_item["name"]
                    source = synonym_item["source"]
                    unversioned_name = ""

                    # check if there's any addtional mapping for the source, remap if so
                    if additional_sources_mapping:
                       source = additional_sources_mapping.get(source, source) # use the same name if no matches in additional_sources_mapping dict

                    # try to get unversioned name if applicable 
                    if source in unversionable_sources:
                        unversioned_name = re.sub(r"\.\d+$", "", synonym_name)

                    # put trios if names are not already seen in db 
                    if synonym_name not in synonym_in_db:
                        external_db_id = self.id_from_map_or_die(source, external_db_map, "external_db_map")
                        synonyms_from_json.append( (seq_region_id, synonym_name, external_db_id) )

                    #   put additional unversioned synonyms if there's a sane one
                    if unversioned_name \
                      and unversioned_name != synonym_name \
                      and unversioned_name not in synonym_in_db:
                        synonyms_from_json.append( (seq_region_id, synonym_name, ensembl_internal_synonym_ext_db_id) )

        # run insertion SQL
        self.insert_to_db(
            synonyms_from_json,
            "seq_region_synonym",
            ["seq_region_id", "synonym", "external_db_id"],
            pj(work_dir, "new_seq_region_synonyms"),
            ignore = True
        )


    def add_sr_attribs(seq_region_file: str,
            seq_region_map: dict,
            attrib_type_map: dict,
            work_dir,
            unversion: bool = unversion):
        """
        Add seq_region_attrib(s) from the seq_region_file meta data file. Explicit list is taken from "sr_attrib_types" module param.

        Add seq_region_attrib(s) from the schema/seq_region_schema.json compatible meta data file.
        Explicit list is taken from "sr_attrib_types" module param.
        "sr_attrib_types" defines { json_property -> attrib_type.name } map. If the value is dict,
           its keys are treated as "/"-delimetered "json_path" (i.e. "added_sequence/assembly_provider/name").
        No arrays can be processed. Only simple or "flattable" types.
        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file: return

        # technical / optimization. get atttib_type_id(s)
        # create a smaller map with attrib_type_id(s) as values
        properties_to_use = [] # [frozen]set with the top-level "seq_region" properties, that should be processed
        path_attrib_id_map = dict() # { "flatterned/json/paths" : attrib_id_map ))}
        # fill set and map
        for prop, attrib_type in self.param('sr_attrib_types').items():
            # adding high level properties to process
            properties_to_use.append(prop)
            # adding json paths (or properties themselves) to atrrib_type_id map
            if isinstance(attrib_type, dict): # if using json paths (delimeterd with "/")
                for path, inner_attrib_type in attrib_type.items():
                    path_attrib_id_map[path] = self.id_from_map_or_die(inner_attrib_type, attrib_type_map, "attrib_type_map")
            else:
                path_attrib_id_map[prop] = self.id_from_map_or_die(attrib_type, attrib_type_map, "attrib_type_map")
        # return if there's nothing to add
        if not properties_to_use: return

        properties_to_use = frozenset(properties_to_use)

        # load attributes from seq_region file
        attrib_trios = [] # [ (seq_region_id, attrib_id, value)... ] list of trios for inserting into db 
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in seq_regions:
                # get seq_region_id (perhaps, by using unversioned name)
                seq_region_name, seq_region_id, unversioned_name = \
                    self.name_and_id_from_seq_region_item(seq_region, try_unversion = unversion)

                # iterate through properties
                for prop_name in propeties_to_use_map:
                    if prop_name not in seq_region:
                        continue
                    # flattern path
                    path_attrib_id_values_list = self.flattern_seq_region_item(seq_region, prop_name, path_attrib_id_map)
                    # fill attrib_trios
                    for (path, attrib_id, value) in path_attrib_id_values_list:
                        attrib_trios.append( (seq_region_id, attrib_id, value) )

        # run insertion SQL
        self.insert_to_db(
            attrib_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            pj(work_dir, "brc4_ebi_seq_region_synonyms"),
            ignore = True
        )


    def flattern_seq_region_item(self, seq_region: dict, prop_name: str; path_attrib_id_map: dict, sep: str = "/") -> list:
        """
        Flattern seq_region[property] and store corresponding [ (json_path, attrib_id, value)... ] (as list of trios).

        Only works for simple properties or dicts with no arrays on the path. Basically, implemets tree traversal.
        Utility function used by the `add_sr_attribs` method
        """
        res = []
        # is there anything to do
        if prop_name not in seq_region: return res

        # set up
        value = seq_region[prop_name]
        paths_to_go = [ (prop_name, value) ] # storing path and the corresponding value, to prevent repetetive traversals
        # iterate
        while paths_to_go:
            (path, value) = paths_to_go.pop() # get last item
            if isinstance(value, list) :
               # perhaps, it's better to raise exception then to continue silently
               continue
            if isinstance(value, dict) :
                # if value is a complex object, add its leaves
                for key, val in value.items():
                    paths_to_go.append( (f"{path}{sep}{key}", val) )
                continue
            # if value is simple
            attrib_id = path_attrib_id_map.get(path, None)
            if attrib_id:
                res.append( (path, attrib_id, value) )
        # return what ever we have
        return res


    def add_sr_ebi_brc4_names(self,
                              seq_region_file: str,
                              seq_region_map: dict,
                              attrib_type_map: dict,
                              work_dir: str,
                              unversion: bool = False):
        """
        Add "(EBI|BRC4)_seq_region_name" seq_region_attrib(s) either from the seq_region_file meta data file, or from original seq_egion names.

        Add "(EBI|BRC4)_seq_region_name" seq_region_synonym from the schema/seq_region_schema.json compatible meta data file or from the original seq_region_names.
        A special case of attributes adding with default values derived from seq_region names.
        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file: return

        # technical / optimization. get atttib_type_id(s) for "(EBI|BRC4)_seq_region_name"
        tagged_sr_name_attrib_id = {
            tag : self.id_from_map_or_die(f"{tag}_seq_region_name", attrib_type_map, "attrib_type_map") for tag in ["EBI", "BRC4"]
        }

        # load BRC4/EBI name from seq_region file
        brc4_ebi_name_attrib_trios = [] # [ (seq_region_id, attrib_id, value)... ] list of trios for inserting into db 
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in seq_regions:
                # get seq_region_id (perhaps, by using unversioned name)
                seq_region_name, seq_region_id, unversioned_name = \
                    self.name_and_id_from_seq_region_item(seq_region, try_unversion = unversion)
                # append attribs to the brc4_ebi_name_attrib_trios list
                for tag in ["BRC4", "EBI"]:
                    attrib_name = f"{tag}_seq_region_name"
                    attrib_id = tagged_sr_name_attrib_id[tag]
                    value = seq_region.get(attrib_name, seq_region_name)
                    brc4_ebi_name_attrib_trios.append( (seq_region_id, attrib_id, value) )

        # run insertion SQL
        self.insert_to_db(
            brc4_ebi_name_attrib_trios,
            "seq_region_attrib",
            ["seq_region_id", "attrib_type_id", "value"],
            pj(work_dir, "brc4_ebi_seq_region_synonyms"),
            ignore = True
        )


    def add_karyotype_data(self,
                           seq_region_file: str,
                           seq_region_map: dict,
                           attrib_type_map: dict,
                           work_dir: str,
                           unversion: bool = False):
        """
        Adds various karyotypic data from seq_region file and assembly metadata (if present).

        Adds various karyotypic data from the schema/seq_region_schema.json compatible meta data file and assembly metadata (if present).

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible
        """
        # add karyotyope bands data
        regions_with_karyotype_bands = self.add_kayryotype_bands(seq_region_file,
                                                                 seq_region_map,
                                                                 attrib_type_map,
                                                                 pj(work_dir, "karyotype_bands"),
                                                                 unversion = unversion)

        # add karyotype ranks attributes
        # HERE
        asm_meta = self.from_param("genome_data","assembly")
        add_cs_tag = self.param("cs_tag_for_ordered")
        self.add_chr_karyotype_rank(asm_meta, pj(wd,"karyotype"), add_cs_tag)


    def add_karyotype_bands(self,
                            seq_region_file: str,
                            seq_region_map: dict,
                            attrib_type_map: dict,
                            work_dir: str,
                            unversion: bool = False,
                            karyotype_bands_property = "karyotype_bands") -> list:
        """
        Add karyotypic data from the seq_region metafile.

        Add karyotypic data from the schema/seq_region_schema.json compatible meta data file.
        Returns list of [ (seq_region_name, seq_region_id, unversioned_name) ] trios for seq_regions having karyotype bands info.

        If unversion is true:
          * the unversioned synonym would be used to get the seq_region_id from "seq_region_map" if possible

        Too close to the DB schema.
        """
        os.makedirs(work_dir, exist_ok=True)

        # return if there's nothing to add
        if not seq_region_file: return

        # load BRC4/EBI name from seq_region file
        band_tuples = [] # [ (seq_region_id, seq_region_start, seq_region_end, band|"NULL", stain|"NULL")... ] list of tuples for inserting into db 
        with open(seq_region_file) as in_file:
            seq_regions = list(json.load(in_file))
            for seq_region in filter(lambda sr: karyotype_bands_property in sr, seq_regions):
                # iterate through all seq_regions having "karyotype_bands" 

                # get seq_region_id (perhaps, by using unversioned name)
                seq_region_name, seq_region_id, unversioned_name = \
                    self.name_and_id_from_seq_region_item(seq_region, try_unversion = unversion)

                # append bands to the band_tuples list
                for band in seq_region[ karyotype_bands_property ]:
                    # print("BAND: " + str(band), file = sys.stderr)
                    # coords
                    seq_region_start = band["start"]
                    seq_region_end = band["end"]
                    # band_name and stain
                    band_name = self.quote_or_null( band.get("name", None) )
                    stain = self.quote_or_null( band.get("stain", None) )
                    # special cases for stain
                    structure = band.get("structure", None)
                    if structure == "telomere":
                        stain = self.quote_or_null("TEL")
                    elif structure == "centromere":
                        stain = self.quote_or_null("ACEN")

                    # append tuple
                    band_tuples.append( (seq_region_id, eq_region_start, seq_region_end, band_name, stain) )

        # run insertion SQL
        self.insert_to_db(
            band_tuples,
            "karyotype",
            ["seq_region_id", "seq_region_start", "seq_region_end", "band", "stain"]
            pj(work_dir, "karyotype_insertion"),
            ignore = True
        )


    def quote_or_null(self, val: str, quotes: str = "'", null: str = "NULL") -> str;
        """
        Return `val` wrapped in `quotes` or `null` value
        """
        if val is None: return null
        return f"{quotes}{val}{quotes}"


    # STAGES
    def add_chr_karyotype_rank(self, meta, wd, add_cs_tag = None):
        # get order from  meta["chromosome_display_order"] , omit unmentioned
        #   otherwise get toplevel "chromosome" seq_regions, sort by seq_region_id
        os.makedirs(wd, exist_ok=True)
        sr_ids = []
        tag = "chromosome_display_order"
        chr_order = meta and tag in meta and meta[tag] or None
        if (chr_order == None):
            # get chromosome id and karyotype_rank id
            ids_sql = r'''select distinct sr.seq_region_id as seq_region_id
                        from seq_region sr, seq_region_attrib sra, coord_system cs, attrib_type at
                        where sr.seq_region_id = sra.seq_region_id
                          and sr.coord_system_id = cs.coord_system_id
                          and sra.attrib_type_id = at.attrib_type_id
                          and (
                               ( cs.name = "chromosome" and at.code = "toplevel" )
                            or ( at.code = "coord_system_tag" and sra.value = "chromosome" )
                          )
                       order by seq_region_id
                      ;'''
            ids_log_pfx = pj(wd,'chr_ids')
            self.run_sql_req(ids_sql, ids_log_pfx)
            # load
            with open(ids_log_pfx + ".stdout") as f:
                for line in f:
                    if (line.startswith("seq_region_id")):
                        continue
                    (sr_id, ) = line.strip().split("\t")
                    sr_ids.append(sr_id)
            sr_ids = [ (_id, i) for i, _id in enumerate(sr_ids, start = 1) ]
        else:
            # show only chromosomes from the chromosome_order
            # get names, syns from db
            chr_rank = { name : rank for rank, name in enumerate(chr_order, start = 1) }
            # get syns
            syns_out_pfx = pj(wd, "syns_from_core")
            self.get_db_syns(syns_out_pfx)
            # load into dict
            with open(syns_out_pfx + ".stdout") as syns_file:
                for line in syns_file:
                    (sr_id, name, syn) = line.strip().split("\t")
                    for _name in [name, syn]:
                        if _name in chr_rank:
                           sr_ids.append((int(sr_id), chr_rank[_name]))
            sr_ids=list(set(sr_ids))
            # assert chr_rank is not reused
            if len(sr_ids) != len(frozenset(map(lambda p: p[1], sr_ids))):
                raise Exception("karyotype_rank is reused: %s" % (str(sr_ids)))
            # assert seq_region_id is not reused
            if len(sr_ids) != len(frozenset(map(lambda p: p[0], sr_ids))):
                raise Exception("same seq_region with different karyotype_rank or wrong seq_regions used in \"chromosome_display_order\". known: %s" % (str(sr_ids)))
            # trying to set chromosome tag
            #   should not change or add if seq_region_tag is already loaded (INSERT IGNORE used)
            if len(sr_ids) > 0 and add_cs_tag is not None:
              tag = "coord_system_tag"
              sr_ids_chr = [ (_id, add_cs_tag) for _id, _  in sr_ids ]
              self.set_sr_attrib(tag, sr_ids_chr, pj(wd, "sr_attr_set_"+tag))

        # insert attrib sql
        if len(sr_ids) > 0:
            tag = "karyotype_rank"
            self.set_sr_attrib(tag, sr_ids, pj(wd, "sr_attr_set_"+tag))


    def unversion_scaffolds(self, cs_rank, logs):
        """
        Unversion scaffold, remove ".\d$" from seq_region.names if there's a need

        Non-versioned syns for contigs (lower, sequence level), versioned for the rest.
        """
        seq_cs, max_rank = max([ (c, r) for c, r in cs_rank.items()], key = lambda k: k[1])
        for cs in cs_rank:
            if cs == seq_cs:
                # for non-sequence level cs, store original name (with version) as "sr_syn_src" synonym
                xdb = self.param("sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, pj(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region_synonym", "synonym", pj(logs, "unv_srs", cs))
            else:
                # for sequence-level cs, store original name (with version) as "versioned_sr_syn_src" synonym
                xdb = self.param("versioned_sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, pj(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region", "name", pj(logs, "unv_sr", cs))


    def coord_sys_order(self, cs_order_str):
        cs_order_lst = map(lambda x: x.strip(), cs_order_str.split(","))
        return { e:i for i,e in enumerate(filter(lambda x: len(x)>0, cs_order_lst)) }


    def used_cs_ranks(self, agps, cs_order, noagp_default = None):
        # rank cs_names, met in agps.keys ("-" separated), i.e. "scaffold-contig"
        #   only agps keys used, values are ignored
        #   use noagp_cs_name_default for "noagp" assemblies
        if agps is None:
            if noagp_default is None:
                raise Exception("NoAGP assembly with no default coordinate system name")
            cs_used_set = frozenset([noagp_default])
        else:
            cs_used_set = frozenset(sum(map(lambda x: x.split("-"), agps.keys()),[]))

        cs_unknown = cs_used_set.difference(cs_order.keys())
        if (len(cs_unknown) > 0):
            raise Exception("Unknown coordinate system(s) %s" % {str(cs_unknown)})
        return { e:i for i,e in enumerate(sorted(cs_used_set,key=lambda x:-cs_order[x]), start=1) }


    def prune_agps(self, agps, cs_order, agps_pruned_dir, pruning = True):
        # when loading agp sort by:
        #   highest component (cmp) cs level (lowest rank)
        #   lowest difference between cs ranks (asm - cmp)
        #   i.e: chromosome-scaffold scaffold-chunk chromosome-chunk
        #   if no agps return empty pruned result
        if agps is None: return None

        agp_cs_pairs = list(map(lambda x: [x]+x.split("-"), agps.keys()))

        agp_levels = [ (x[0], cs_order[x[1]], cs_order[x[2]]) for x in agp_cs_pairs ]
        bad_agps = list(filter(lambda x: x[1] < x[2], agp_levels))
        if (len(bad_agps) > 0):
            raise Exception("component cs has higher order than assembled cs %s" % (str(bad_agps)))

        agp_levels_sorted = [ e[0] for e in sorted(agp_levels, key=lambda x:(-x[2], x[1]-x[2])) ]

        #prune agps
        agps_pruned = {}
        used_components = set()
        if not pruning:
            used_components = None
        for asm_cmp in agp_levels_sorted:
            agp_file_src = agps[asm_cmp]
            agp_file_dst = pj(agps_pruned_dir, asm_cmp + ".agp")
            if self.agp_prune(agp_file_src, agp_file_dst, used_components) > 0:
                agps_pruned[asm_cmp] = agp_file_dst
        return agps_pruned


    def load_seq_data(self, fasta, agps, cs_rank, log_pfx):
        """loads sequence data for various coordinate systems accordingly with their rank"""
        asm_v = self.asm_name()

        sequence_rank = max(cs_rank.values())
        for (cs, rank) in sorted(cs_rank.items(), key=lambda p: -p[1]):
           logs = pj(log_pfx, "%02d_%s" %(rank, cs) )
           if (rank == sequence_rank):
               self.load_cs_data(cs, rank, "fasta", asm_v, fasta, logs, loaded_regions = None, seq_level = True)
           else:
               useful_agps = list(filter(lambda x: cs in x, agps and agps.keys() or []))
               if len(useful_agps) == 0:
                   raise Exception("non-seq_level cs %s has no agps to assemble it from" % (cs))
               loaded_regions = set()
               for pair, agp_file_pruned in map(lambda k: (k, agps[k]), useful_agps):
                   if (not pair.startswith(cs+"-")):
                       continue
                   self.load_cs_data(cs, rank, pair, asm_v, agp_file_pruned, logs, loaded_regions)

    def load_cs_data(self,
                     cs, rank, pair, asm_v,
                     src_file, log_pfx,
                     loaded_regions = None, seq_level = False):
        """creates a coord_system and loads sequence or assembly(AGP) data for corresponding seqregions

           doesn't load already seen sequences
        """
        # NB load_seq_region.pl and load_agp.pl are not failing on parameter errors (0 exit code)
        os.makedirs(dirname(log_pfx), exist_ok=True)
        if seq_level:
            self.load_seq_region(cs, rank, asm_v, src_file, log_pfx, seq_level)
        elif loaded_regions is not None:
            new_regions = set()
            clean_file = src_file + ".regions_deduped"
            self.filter_already_loaded_regions_from_agp(src_file, clean_file, loaded_regions, new_regions)
            self.load_seq_region(cs, rank, asm_v, clean_file, log_pfx, seq_level)
            loaded_regions.update(new_regions)
        if not seq_level:
            self.load_agp(pair, asm_v, src_file, log_pfx)


    def filter_already_loaded_regions_from_agp(self, src_file, dst_file, loaded_regions, new_regions):
       with open(src_file) as src:
           with open(dst_file, "w") as dst:
                for line in src:
                    fields = line.strip().split("\t")
                    ( asm_id, asm_start, asm_end, asm_part,
                      type_,
                      cmp_id, cmp_start, cmp_end, cmp_strand
                    ) = fields
                    if type_ in "NU" or asm_id in loaded_regions:
                        continue
                    new_regions.add(asm_id)
                    print(line.strip(), file = dst)


    def agp_prune(self, from_file: str, to_file: str, used: set = None):
        """
        Remove already components from the AGP file if they are seen in "used" set
        """
        # reomve used component
        #   and GAPS as they are not used by 'ensembl-analysis/scripts/assembly_loading/load_agp.pl'
        os.makedirs(dirname(to_file), exist_ok=True)
        open_ = self.is_gz(from_file) and gzip.open or open
        if used is None:
            cmd = r'''{_cat} {_file} > {_out}'''.format(
                _cat = self.is_gz(from_file) and "zcat" or "cat",
                _file = from_file,
                _out = to_file
            )
            print("running %s" % (cmd), file = sys.stderr)
            sp.run(cmd, shell=True, check=True)
            return 1
        writes = 0
        with open_(from_file, "r") as src:
            with open(to_file, "w") as dst:
                for line in src:
                    fields = line.strip().split("\t")
                    ( asm_id, asm_start, asm_end, asm_part,
                      type_,
                      cmp_id, cmp_start, cmp_end, cmp_strand
                    ) = fields
                    if type_ in "NU" or cmp_id in used:
                        continue
                    used.add(cmp_id)
                    print(line.strip(), file = dst)
                    writes += 1
        return writes


    def get_external_db_mapping(self) -> dict:
        """
        Get a map from a file for external_dbs to Ensembl dbnames from "external_db_map" module(!) param
        """
        external_map_path = self.param("external_db_map")
        db_map = dict()
        if external_map_path is None: return db_map

        # Load the map
        with open(external_map_path, "r") as map_file:
            for line in map_file:
                if line.startswith("#"): continue
                line = re.sub(r'#.*', '', line)
                if re.match(r'^\s*$', line): continue
                (from_name, to_name, *rest) = line.strip().split("\t")
                if len(rest) > 0 and rest[0].upper() != "SEQ_REGION": continue
                if to_name == "_IGNORE_": continue
                db_map[from_name] = to_name
        return db_map


    # UTILS
    def db_string(self):
        return "-dbhost {host_} -dbport {port_} -dbuser {user_} -dbpass {pass_} -dbname {dbname_} ".format(
            host_ = self.param("dbsrv_host"),
            port_ = self.param("dbsrv_port"),
            user_ = self.param("dbsrv_user"),
            pass_ = self.param("dbsrv_pass"),
            dbname_ = self.param("db_name")
        )


    def is_gz(self, filename):
      return filename.endswith(".gz")


    def asm_name(self):
        asm = self.from_param("genome_data","assembly")
        if "name" not in asm:
            raise Exception("no assembly/name in genome_data")
        return asm["name"]


    # TODO: add some metafunc setter getter
    def from_param(self, param, key, not_throw = False):
        data = self.param_required(param)
        if key not in data:
            if not_throw:
                return None
            else:
                raise Exception("Missing required %s data: %s" % (param , key))
        return data[key]


    def param_bool(self, param):
        val = self.param(param)
        return bool(val) and "0" != val


    def load_map_from_sql_stdout(self, in_file, skip_header = False):
        """
        Load map from the SQL output

        Process input in_file with "key  value" pairs
          and load then into the {key : value} map.
        Skips header if skip_heade.
        """
        with open(in_file) as pairs_file:
            for line in pairs_file:
                if skip_header:
                    skip_header = False
                    continue
                (key, val) = line.strip().split("\t")
                data[key] = val
        return data


    def name_and_id_from_seq_region_item(self, seq_region_item: dict, seq_region_map: dict, try_unversion: bool = False, throw_missing: bool = True) -> (str, str, str):
        """
        Get (seq_region_name, seq_region_id, unversioned_name) from seq_region_item struct(dict)

        Gets unversioned_name only if "try_unversion" is True.
        Throws exception if not able to get seq_region_id from "seq_region_map" and "throw_missing" is true.
        """
        #   get seq_region_id (perhaps, by using unversioned name)
        seq_region_name = seq_region["name"]
        seq_region_id = seq_region_map.get(seq_region_name, None)
        if seq_region_id is None and unversion:
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
        cmd = r'''{_cat} {_file} | sed -r '/^[^>]/ {{ s/[{_IUPAC}]/N/g; s/{_iupac}/n/g }}' > {_out}'''.format(
            _cat = self.is_gz(from_file) and "zcat" or "cat",
            _file = from_file,
            _IUPAC = IUPAC.upper(),
            _iupac = IUPAC.lower(),
            _out = to_file
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)


    def load_seq_region(self, cs: str, rank: str, asm_v: str, src_file: str, log_pfx: str, seq_level = False):
        """ensembl-analysis script (load_seq_region.pl) based utility for loading seq_regions FASTA sequences"""
        en_root = self.param_required("ensembl_root_dir")
        cmd = (r'''{_loader} {_db_string} -coord_system_version {_asm_v} -default_version ''' +
               r'''    -rank {_rank} -coord_system_name {_cs} {_sl_flag} -{_tag}_file {_file}''' +
               r'''     > {_log}.stdout 2> {_log}.stderr''').format(
            _loader = "perl %s" % (pj(en_root, r"ensembl-analysis/scripts/assembly_loading/load_seq_region.pl")),
            _db_string = self.db_string(),
            _asm_v = asm_v,
            _rank = rank,
            _cs = cs,
            _sl_flag = seq_level and "-sequence_level" or "",
            _tag = seq_level and "fasta" or "agp",
            _file = src_file,
            _log = "%s_seq" % (log_pfx),
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)


    def load_agp(self, pair, asm_v, src_file, log_pfx):
        """ensembl script (load_agp.pl) based utility for loading seq_regions assembly data (AGPs)"""
        en_root = self.param_required("ensembl_root_dir")
        (asm_n, cmp_n) = pair.strip().split("-")
        cmd = (r'''{_loader} {_db_string} -assembled_version {_asm_v} ''' +
               r'''    -assembled_name {_asm} -component_name {_cmp} ''' +
               r'''    -agp_file {_file} ''' +
               r'''    > {_log}.stdout 2> {_log}.stderr''').format(
            _loader = "perl %s" % (pj(en_root, r"ensembl-analysis/scripts/assembly_loading/load_agp.pl")),
            _db_string = self.db_string(),
            _asm_v = asm_v,
            _asm = asm_n,
            _cmp = cmp_n,
            _file = src_file,
            _log = "%s_agp_%s" % (log_pfx, pair.replace("-","_")),
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)


    def set_toplevel(self, log_pfx, ignored_cs = []):
        """
        Set toplevel(6) seq_region_attrib using ensembl script.

        Uses set_toplevel.pl ensembl script.
        """
        # set top_level(6) seq_region_attrib
        os.makedirs(dirname(log_pfx), exist_ok=True)
        en_root = self.param_required("ensembl_root_dir")
        cmd = (r'''{_set_tl} {_db_string} {_ignored_cs} ''' +
               r'''     > {_log}.stdout 2> {_log}.stderr''').format(
            _set_tl = "perl %s" % (pj(en_root, r"ensembl-analysis/scripts/assembly_loading/set_toplevel.pl")),
            _db_string = self.db_string(),
            _ignored_cs = " ".join(map(lambda x: "-ignore_coord_system %s" % (x), ignored_cs)),
            _log = log_pfx,
        )
        print("running %s" % (cmd), file = sys.stderr)
        sp.run(cmd, shell=True, check=True)

        # remove toplevel attribute for seq_regions that are components 
        self.remove_components_from_toplevel(log_pfx)


    ## SQL executor and utilities using plain SQL
    def run_sql_req(self, sql, log_pfx, from_file = False):
        os.makedirs(dirname(log_pfx), exist_ok=True)
        en_root = self.param_required("ensembl_root_dir")

        sql_option = r''' -sql '{_sql}' '''.format(_sql = sql)
        if from_file:
            sql_option = r''' < '{_sql}' '''.format(_sql = sql)

        cmd = r'''{_dbcmd} -url "{_srv}{_dbname}" {_sql_option} > {_out} 2> {_err}'''.format(
            _dbcmd = 'perl %s/ensembl-hive/scripts/db_cmd.pl' %(en_root),
            _srv = self.param("dbsrv_url"),
            _dbname = self.param("db_name"),
            _sql_option = sql_option,
            _out = log_pfx + ".stdout",
            _err = log_pfx + ".stderr"
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)


    def add_contig_ena_attrib(self, log_pfx, cs_name = "contig"):
        """
        Add ENA attrib for contigs if their names are ENA accessions

        Nno sequence_level checks are used -- just cs name.
        See ensembl-datacheck/lib/Bio/EnsEMBL/DataCheck/Checks/SeqRegionNamesINSDC.pm .
        SQL code.
        """
        sql = r'''insert into seq_region_attrib (seq_region_id, attrib_type_id, value)
                select
                  sr.seq_region_id, at.attrib_type_id, "ENA"
                from
                  seq_region sr, coord_system cs, attrib_type at
                where   sr.coord_system_id = cs.coord_system_id
                    and cs.name = "%s"
                    and at.code = "external_db"
              ;''' % (cs_name)
        return self.run_sql_req(sql, log_pfx)


    def copy_sr_name_to_syn(self, cs, x_db, log_pfx):
        """
        Store original seq_region names as seq_region_synonym

        Store original seq_region names (from a given cood_systen, "cs" param) as seq_region_synonyms (using "x_db" external source name)
        SQL code.
        """
        asm_v = self.asm_name()
        sql = r'''insert into seq_region_synonym (seq_region_id, synonym, external_db_id)
                  select
                      sr.seq_region_id, sr.name, xdb.external_db_id
                  from
                     seq_region sr, external_db xdb, coord_system cs
                  where   xdb.db_name = "%s"
                      and sr.coord_system_id = cs.coord_system_id
                      and cs.name = "%s"
                      and cs.version = "%s"
                      and sr.name like "%%._"
                ;''' % (x_db, cs, asm_v)
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
            sql = r'''insert ignore into meta (species_id, meta_key, meta_value) values
                    (1, "assembly.mapping", "{_higher}:{_v}|{_lower}:{_v}")
                  ;'''.format(_v = asm_v, _higher = higher, _lower = lower)
            self.run_sql_req(sql, pj(log_pfx, pair))

    def remove_components_from_toplevel(self, log_pfx):
        """
        Remove toplevel attribute for seq_regions that are "components" (parts of different seq_regions).

        SQL code.
        """

        # get list of seq_regions that are components
        sql_not_toplevel_list = r'''
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
        '''
        # perhaps, make sense to check sr_c.coord_system_id != sr_a.coord_system_id
        self.run_sql_req(sql_not_toplevel_list, ".".join([log_pfx, "not_toplevel_list"]))

        # delete wrongly assigned attribs
        sql_not_toplevel_delete = r'''
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
        '''
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
        sql = r'''update {_tbl} t, seq_region sr, coord_system cs
                    set
                      t.{_fld} = substr(t.{_fld},  1, locate(".", t.{_fld}, length(t.{_fld})-2)-1)
                    where t.{_fld} like "%._"
                      and t.seq_region_id = sr.seq_region_id
                      and sr.coord_system_id = cs.coord_system_id
                      and cs.name = "{_cs}"
                      and cs.version = "{_asm_v}"
                ;'''.format(
                    _tbl = tbl,
                    _fld = fld,
                    _cs = cs,
                    _asm_v = asm_v
                )
        return self.run_sql_req(sql, log_pfx)


    def nullify_ctg_cs_version(self, log_pfx: str):
        """
        Nullify every CS version with rank larger than that of "contig", but don't nullify toplevel ones.

        SQL code
        """
        asm_v = self.asm_name()
        # get cs_info (and if they have toplevel regions)
        sql = r'''select cs.coord_system_id as coord_system_id,
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
              ;'''.format(_asm_v = asm_v)
        # run_sql
        toplvl_pfx = pj(log_pfx,"toplvl_info")
        self.run_sql_req(sql, toplvl_pfx)
        # load info
        cs_info = []
        with open(toplvl_pfx + ".stdout") as f:
            header = next(f).strip().split("\t")
            for line in f:
                cs_info.append(dict(zip(header, line.strip().split())))
        # choose cs rank threshold to start clearing version from
        seq_rank = max(map(lambda cs: int(cs["rank"]), cs_info))
        nullify_cs_version_from = self.param("nullify_cs_version_from")
        ctg_lst = list(filter(lambda cs: cs["name"] == nullify_cs_version_from, cs_info))
        clear_thr = ctg_lst and int(ctg_lst[0]["rank"]) or seq_rank
        clear_lst = [ (cs["coord_system_id"], cs["name"]) for cs in cs_info
                        if (bool(int(cs["no_toplevel"])) and int(cs["rank"]) >= clear_thr) ]
        # run sql
        if clear_lst:
            clear_pfx = pj(log_pfx, "clear")
            with open(clear_pfx + ".sql", "w") as clear_sql:
                for (cs_id, cs_name) in clear_lst:
                    sql = r'''
                        update meta set
                            meta_value=replace(meta_value, "|{_cs_name}:{_asm_v}", "|{_cs_name}")
                            where meta_key="assembly.mapping";
                        update coord_system set version = NULL where coord_system_id = {_cs_id};
                    '''.format(_asm_v = asm_v, _cs_name = cs_name, _cs_id = cs_id)
                    print(sql, file = clear_sql)
            self.run_sql_req(clear_pfx + ".sql", clear_pfx, from_file = True)


    def load_map_from_core_db(self, table, cols, work_dir) -> dict:
        """
        Load 2 "cols" from core db "table" as map

        Load { cols[0] : cols[1] } map from the core db "table"
        SQL code
        """
        out_pfx = pj(work_dir, f"{table}_map")
        sql = f'''select {cols[0]}, {cols[1]} FROM {table};'''
        res = self.run_sql_req(sql, out_pfx)

        out_file = out_pfx + ".stdout"
        data = self.load_map_from_sql_stdout(out_file, skip_header = True)
        if not data:
            raise Exception(f"No '{table}' map loaded from '{out_file}'")
        return data


    def load_seq_region_synonyms_trios_from_core_db(self, work_dir: str) -> list:
        # was get_db_syns
        """
        Load seq_region_synonyms from from core db into [(seq_region_id, name, synonym)...] list

        SQL code
        """
        out_pfx = pj(work_dir, f"seq_region_synonyms")
        sql = r'''select sr.seq_region_id as seq_region_id, sr.name, srs.synonym
                 from seq_region sr left join seq_region_synonym srs
                 on sr.seq_region_id = srs.seq_region_id
                 order by sr.seq_region_id
              ;'''
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
            self,
            list_of_tuples: list,
            table_name: str,
            col_names: list,
            work_dir: str,
            ignore: bool = True
        )
        """
        Insert into the core db's {table_name} tuples from {list_of_tuples} as col_names.

        SQL code
        """
        # return if nothing to do
        if not list_of_tuples: return

        # prepare request parts
        ignore_str = ignore and "IGNORE" or ""
        cols_str = ", ".join(col_names)

        # generate file with the insert SQL command
        insert_sql_file = pj(work_dir, "insert.sql")
        with open(insert_sql_file, "w") as sql:
            print("INSERT {ignore_str} INTO {table_name} ({col_names}) VALUES", file=sql)
            values_sep = ""
            for tpl in list_of_tuple:
                tpl_str = ", ".join(map v: str(v), tpl)
                print(f"{values_sep}({tpl_str})", file = sql)
                values_sep = ", "
            print(";", file=sql)

        # run insert SQL from file
        self.run_sql_req(insert_sql_file, pj(wd, "insert_syns"), from_file = True)


