#!env python3

from BCBio import GFF

import eHive
import gzip
import json
import os
import subprocess as sp
import sys

from os.path import dirname, join as pj

from Bio import SeqIO

import io
from math import floor

class LoadSequenceData(eHive.BaseRunnable):

    def param_defaults(self):
        return {
            'cs_order' : 'chunk,contig,supercontig,non_ref_scaffold,scaffold,superscaffold,linkage_group,chromosome',
            'IUPAC' : 'RYKMSWBDHV',
            'unversion_scaffolds' : 0,
            'versioned_sr_syn_src' : 'INSDC', # 50710
            'sr_syn_src' : 'VB_Community_Symbol', # 211
            'sr_attrib_types' : {
                'circular' : 'circular_seq',
                'codon_table' : 'codon_table',
                'location' : 'SO_term',
                'non_ref' : 'non_ref',
                'karyotype_bands' : 'karyotype_bands',
            },
            'not_toplevel_cs' : [], # i.e. "contig", "non_ref_scaffold"
            'nullify_cs_version_from' : 'contig',
        }


    def run(self):
        # params
        en_root = self.param_required("ensembl_root_dir")
        wd = self.param_required("work_dir")

        # initialize whatever
        genome = self.from_param("manifest_data", "genome")
        sra = self.from_param("manifest_data", "seq_region")

        # TODO
        # split into contigs, add AGP
        # load data with no agps ??? m.b. create empty cs-cs agps
        # omit, split

        # rename IUPAC to N symbols using sed
        fasta_raw = self.from_param("manifest_data", "fasta_dna")
        fasta_clean = pj(wd, "fasta", "seq_no_iupac.fasta")
        self.remove_IUPAC(fasta_raw, fasta_clean)

        agps = self.from_param("manifest_data", "agp")
        cs_order = self.coord_sys_order(self.param("cs_order"))
        cs_rank = self.used_cs_ranks(agps, cs_order)

        agps_pruned_dir = pj(wd, "agps_pruned")
        agps_pruned = self.prune_agps(agps, cs_order, agps_pruned_dir, self.param_bool("prune_agp"))

        self.load_seq_data(fasta_clean, agps_pruned, cs_rank, pj(wd, "load"))

        self.add_contig_ena_attrib(pj(wd, "load", "set_ena"))

        unversion_scaffolds = self.param_bool("unversion_scaffolds")
        if unversion_scaffolds:
            self.unversion_scaffolds(cs_rank, pj(wd, "unversion_scaffolds"))

        seq_reg_file = self.from_param("manifest_data", "seq_region")
        self.add_sr_synonyms(seq_reg_file, pj(wd, "seq_region_syns"), unversion_scaffolds)

        self.add_sr_attribs(seq_reg_file, pj(wd, "seq_region_attr"), karyotype_info_tag = "karyotype_bands")

        self.add_asm_mappings(agps_pruned.keys(), pj(wd, "asm_mappings"))

        self.set_toplevel(pj(wd, "set_toplevel"), self.param("not_toplevel_cs"))

        asm_meta = self.from_param("genome_data","assembly")
        self.add_chr_karyotype_rank(asm_meta, pj(wd,"karyotype"))


    # STAGES
    def add_asm_mappings(self, cs_pairs, log_pfx):
        # nullifies asm_mappings contig versions as well, but don't nullify toplevel
        asm_v = self.from_param("genome_data","assembly")["name"]
        for pair in cs_pairs:
            higher, lower = pair.strip().split("-")
            sql = r'''insert ignore into meta (species_id, meta_key, meta_value) values
                    (1, "assembly.mapping", "{_higher}:{_v}|{_lower}:{_v}")
                  ;'''.format(_v = asm_v, _higher = higher, _lower = lower)
            self.run_sql_req(sql, pj(log_pfx, pair))
        self.nullify_ctg_cs_version(pj(log_pfx, "nullify_cs_versions"))


    def nullify_ctg_cs_version(self, log_pfx):
        # nullify every cs with rank larger than contig, but don't nullify toplevel ones
        asm_v = self.from_param("genome_data","assembly")["name"]
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
                        if (bool(cs["no_toplevel"]) and int(cs["rank"]) >= clear_thr) ]
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


    def add_chr_karyotype_rank(self, meta, wd):
        # get order from  meta["chromosome_display_order"] , omit unmentioned
        #   otherwise get toplevel "chromosome" seq_regions, sort by seq_region_id
        os.makedirs(wd, exist_ok=True)
        sr_ids = []
        tag = "chromosome_display_order"
        chr_order = meta and tag in meta and meta[tag] or None
        if (chr_order == None):
            # get chromosome id and karyotype_rank id
            ids_sql = r'''select sr.seq_region_id as seq_region_id
                        from seq_region sr, seq_region_attrib sra, coord_system cs, attrib_type at
                        where sr.seq_region_id = sra.seq_region_id
                          and sr.coord_system_id = cs.coord_system_id
                          and sra.attrib_type_id = at.attrib_type_id
                          and cs.name = "chromosome"
                          and at.code = "toplevel"
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
            order = dict(map(lambda e: (e[1], e[0]), enumerate(chr_order, start = 1)))
            # get syns
            syns_out_pfx = pj(wd, "syns_from_core")
            self.get_db_syns(syns_out_pfx)
            # load into dict
            with open(syns_out_pfx + ".stdout") as syns_file:
                for line in syns_file:
                    (sr_id, name, syn) = line.strip().split("\t")
                    for _name in [name, syn]:
                        if _name in order:
                           sr_ids.append((int(sr_id), order[_name]))
            # assert order is not reused
            if len(sr_ids) != len(frozenset(map(lambda p: p[1], sr_ids))):
                raise Exception("karyotype_rank is reused: %s" % (str(sr_ids)))
            # assert seq_region_id is not reused
            if len(sr_ids) == len(frozenset(map(lambda p: p[0], sr_ids))):
                raise Exception("same seq_region with different karyotype_rank: %s" % (str(sr_ids)))
        # insert attrib sql
        if len(sr_ids) > 0:
            tag = "karyotype_rank"
            self.set_sr_attrib(tag, sr_ids, pj(wd, "sr_attr_set_"+tag))


    def set_toplevel(self, log_pfx, ignored_cs = []):
        # set top_level(6) seq_region_attrib
        os.makedirs(dirname(log_pfx), exist_ok=True)
        en_root = self.param_required("ensembl_root_dir")
        cmd = (r'''{_set_tl} {_db_string} {_ignored_cs} ''' +
               r'''     > {_log}.stdout 2> {_log}.stderr''').format(
            _set_tl = "perl %s" % (pj(en_root, r"ensembl-pipeline/scripts/set_toplevel.pl")),
            _db_string = self.db_string(),
            _ignored_cs = " ".join(map(lambda x: "-ignore_coord_system %s" % (x), ignored_cs)),
            _log = log_pfx,
        )
        print("running %s" % (cmd), file = sys.stderr)
        sp.run(cmd, shell=True, check=True)
        # remove toplevel attribute for seq_regions that are parts of different seq_regions
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


    def add_sr_attribs(self, meta_file, wd, karyotype_info_tag = None):
        os.makedirs(wd, exist_ok=True)
        # find interesting attribs in meta_file
        attribs_map = self.param("sr_attrib_types")
        interest = frozenset(attribs_map.keys())
        chosen = dict()
        with open(meta_file) as mf:
            data = json.load(mf)
            if not isinstance(data, list):
                data = [ data ]
            for e in data:
                if interest.intersection(e.keys()):
                    chosen[e["name"]] = [e, -1]
        if len(chosen) <= 0:
            return
        # get names, syns from db
        syns_out_pfx = pj(wd, "syns_from_core")
        self.get_db_syns(syns_out_pfx)
        # load into dict
        with open(syns_out_pfx + ".stdout") as syns_file:
            for line in syns_file:
                (sr_id, name, syn) = line.strip().split("\t")
                for _name in [name, syn]:
                    if _name in chosen:
                        sr_id = int(sr_id)
                        if chosen[_name][1] != -1 and chosen[_name][1] != sr_id:
                            raise Exception(
                                "Same name reused by different seq_regions: %d , %d" % (
                                    chosen[_name][1], sr_id
                                )
                            )
                        chosen[_name][1] = sr_id
        # add seq_region_attribs
        for tag, attr_type in attribs_map.items():
            srlist = list(map(
                lambda p:(p[1], p[0][tag]),
                filter(lambda x: tag in x[0], chosen.values())
            ))
            self.set_sr_attrib(
                attr_type,
                srlist,
                pj(wd, "sr_attr_set_"+tag),
                (karyotype_info_tag and tag == karyotype_info_tag)
            )


    def add_karyotype_bands(self, id_val_lst, log_pfx):
        os.makedirs(dirname(log_pfx), exist_ok=True)
        insert_sql_file = log_pfx + "_insert_karyotype_bands.sql"
        if len(id_val_lst) <= 0:
            return
        with open(insert_sql_file, "w") as sql:
            print("insert into karyotype (seq_region_id, seq_region_start, seq_region_end, band, stain) values", file=sql)
            sep = ""
            for _sr_id, _val in id_val_lst:
                for band in _val:
                    # print("BAND: " + str(band), file = sys.stderr)
                    start = band["start"]
                    end = band["end"]
                    name = "name" in band and "'%s'" % (band["name"]) or "NULL"
                    stain = "stain" in band and "'%s'" % (band["stain"]) or "NULL"
                    if "structure" in band:
                        if band["structure"] == "telomere":
                            stain = "'TEL'"
                        elif band["structure"] == "centromere":
                            stain = "'ACEN'"
                    val_str =  ",".join(list(map(lambda x: str(x), [_sr_id, start, end, name, stain])))
                    print ('%s (%s)' % (sep, val_str), file = sql)
                    sep = ","
            print(";", file=sql)
        # run insert sql
        self.run_sql_req(insert_sql_file, log_pfx, from_file = True)


    def set_sr_attrib(self, attr_type, id_val_lst, log_pfx, karyotype_info = False):
        if karyotype_info: 
            return self.add_karyotype_bands(id_val_lst, log_pfx)
        # generaate sql req for loading
        os.makedirs(dirname(log_pfx), exist_ok=True)
        insert_sql_file = log_pfx + "_insert_attribs.sql"
        if len(id_val_lst) <= 0:
            return
        with open(insert_sql_file, "w") as sql:
            print("insert into seq_region_attrib (seq_region_id, attrib_type_id, value) values", file=sql)
            fst = ""
            for _sr_id, _val in id_val_lst:
                if isinstance(_val, bool):
                    _val = int(_val)
                if isinstance(_val, str):
                    _val = '"%s"' % (_val)
                print ('%s (%s, 0, %s)' % (fst, _sr_id, str(_val)), file = sql)
                fst = ","
            print(";", file=sql)
        # run insert sql
        self.run_sql_req(insert_sql_file, log_pfx, from_file = True)
        # update external_db_id
        sql_update_at = r'''update seq_region_attrib sra, attrib_type at
            set sra.attrib_type_id = at.attrib_type_id
            where sra.attrib_type_id = 0
              and at.code = "%s"
        ;''' % (attr_type)
        self.run_sql_req(sql_update_at, log_pfx+"_update_attr_type")

    def get_db_syns(self, out_pfx):
        # get names, syns from db
        sql = r'''select sr.seq_region_id as seq_region_id, sr.name, srs.synonym
                 from seq_region sr left join seq_region_synonym srs
                 on sr.seq_region_id = srs.seq_region_id
                 order by sr.seq_region_id
              ;'''
        return self.run_sql_req(sql, out_pfx)

    def add_sr_synonyms(self, meta_file, wd, unversioned = False):
        os.makedirs(wd, exist_ok=True)
        # get names, syns from db
        syns_out_pfx = pj(wd, "syns_from_core")
        self.get_db_syns(syns_out_pfx)
        # load them into dict
        sr_ids = dict()
        seen_syns_pre = []
        with open(syns_out_pfx + ".stdout") as syns_file:
            for line in syns_file:
                if (line.startswith("seq_region_id")):
                    continue
                (sr_id, name, syn) = line.strip().split("\t")
                sr_ids[name] = int(sr_id)
                seen_syns_pre.append(name)
                if syn != "NULL":
                    seen_syns_pre.append(syn)
        seen_syns = frozenset(seen_syns_pre)
        # load syns from file
        new_syns = dict()
        with open(meta_file) as mf:
            data = json.load(mf)
            if not isinstance(data, list):
                data = [ data ]
            for e in data:
                if "synonyms" not in e:
                   continue
                es = list(filter(lambda s: s not in seen_syns, e["synonyms"]))
                # do we need unversioned syn as well???
                # Nov 2019: no, we need explicit list
                if len(es) <= 0:
                    continue
                en = e["name"]
                eid = en in sr_ids and sr_ids[en] or None
                if unversioned and eid is None:
                    if en[-2] == ".":
                         en = en[:-2]
                         eid = en in sr_ids and sr_ids[en] or None
                if eid is None:
                    raise Exception("Not able to find seq_region for '%s'" % (e["name"]))
                new_syns[eid] = es
        # generaate sql req for loading
        insert_sql_file = pj(wd, "insert_syns.sql")
        if new_syns:
            with open(insert_sql_file, "w") as sql:
                print("insert into seq_region_synonym (seq_region_id, synonym) values", file=sql)
                fst = ""
                for _sr_id, _sr_syn in sum(map(lambda p: [(p[0], s) for s in p[1]], new_syns.items()), [ ]):
                    print ('%s (%s, "%s")' % (fst, _sr_id, _sr_syn), file = sql)
                    fst = ","
                print(";", file=sql)
            # run insert sql
            self.run_sql_req(insert_sql_file, pj(wd, "insert_syns"), from_file = True)
            # update external_db_id
            sql_update_xdb = r'''update seq_region_synonym srs, external_db xdb
                set srs.external_db_id = xdb.external_db_id
                where srs.external_db_id is NULL
                  and xdb.db_name = "%s"
            ;''' % (self.param("sr_syn_src"))
            self.run_sql_req(sql_update_xdb, pj(wd, "update_syns_xdb"))


    def unversion_scaffolds(self, cs_rank, logs):
        # non-versioned syns for contigs, versioned for the rest
        seq_cs, max_rank = max([ (c, r) for c, r in cs_rank.items()], key = lambda k: k[1])
        for cs in cs_rank:
            if cs == seq_cs:
                xdb = self.param("sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, pj(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region_synonym", "synonym", pj(logs, "unv_srs", cs))
            else:
                xdb = self.param("versioned_sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, pj(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region", "name", pj(logs, "unv_sr", cs))


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


    def add_contig_ena_attrib(self, log_pfx):
        # Add ENA attrib for contigs (no sequence_level checks -- just cs name)
        #   (see ensembl-datacheck/lib/Bio/EnsEMBL/DataCheck/Checks/SeqRegionNamesINSDC.pm)
        sql = r'''insert into seq_region_attrib (seq_region_id, attrib_type_id, value)
                select
                  sr.seq_region_id, at.attrib_type_id, "ENA"
                from
                  seq_region sr, coord_system cs, attrib_type at
                where   sr.coord_system_id = cs.coord_system_id
                    and cs.name = "contig"
                    and at.code = "external_db"
              ;'''
        return self.run_sql_req(sql, log_pfx)


    def copy_sr_name_to_syn(self, cs, x_db, log_pfx):
        asm_v = self.from_param("genome_data","assembly")["name"]
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


    def sr_name_unversion(self, cs, tbl, fld, log_pfx):
        # select synonym, substr(synonym,  1, locate(".", synonym, length(synonym)-2)-1)
        #     from seq_region_synonym  where synonym like "%._"
        asm_v = self.from_param("genome_data","assembly")["name"]
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


    def coord_sys_order(self, cs_order_str):
        cs_order_lst = map(lambda x: x.strip(), cs_order_str.split(","))
        return { e:i for i,e in enumerate(filter(lambda x: len(x)>0, cs_order_lst)) }


    def used_cs_ranks(self, agps, cs_order):
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
        asm_v = self.from_param("genome_data","assembly")["name"]

        sequence_rank = max(cs_rank.values())
        for (cs, rank) in sorted(cs_rank.items(), key=lambda p: -p[1]):
           logs = pj(log_pfx, "%02d_%s" %(rank, cs) )
           if (rank == sequence_rank):
               self.load_cs_data(cs, rank, "fasta", asm_v, fasta, logs, loaded_regions = None, seq_level = True)
           else:
               useful_agps = list(filter(lambda x: cs in x, agps.keys()))
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

    def load_seq_region(self, cs, rank, asm_v, src_file, log_pfx, seq_level = False):
        en_root = self.param_required("ensembl_root_dir")
        cmd = (r'''{_loader} {_db_string} -coord_system_version {_asm_v} -default_version ''' +
               r'''    -rank {_rank} -coord_system_name {_cs} {_sl_flag} -{_tag}_file {_file}''' +
               r'''     > {_log}.stdout 2> {_log}.stderr''').format(
            _loader = "perl %s" % (pj(en_root, r"ensembl-pipeline/scripts/load_seq_region.pl")),
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
        en_root = self.param_required("ensembl_root_dir")
        (asm_n, cmp_n) = pair.strip().split("-")
        cmd = (r'''{_loader} {_db_string} -assembled_version {_asm_v} ''' +
               r'''    -assembled_name {_asm} -component_name {_cmp} ''' +
               r'''    -agp_file {_file} ''' +
               r'''    > {_log}.stdout 2> {_log}.stderr''').format(
            _loader = "perl %s" % (pj(en_root, r"ensembl-pipeline/scripts/load_agp.pl")),
            _db_string = self.db_string(),
            _asm_v = asm_v,
            _asm = asm_n,
            _cmp = cmp_n,
            _file = src_file,
            _log = "%s_agp_%s" % (log_pfx, pair.replace("-","_")),
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)


    def db_string(self):
        return "-dbhost {host_} -dbport {port_} -dbuser {user_} -dbpass {pass_} -dbname {dbname_} ".format(
            host_ = self.param("dbsrv_host"),
            port_ = self.param("dbsrv_port"),
            user_ = self.param("dbsrv_user"),
            pass_ = self.param("dbsrv_pass"),
            dbname_ = self.param("db_name")
        )


    def agp_prune(self, from_file, to_file, used = None):
        # reomve used component
        #   and GAPS as they are not used by 'ensembl-pipeline/scripts/load_agp.pl'
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


    def remove_IUPAC(self, from_file, to_file):
        IUPAC = self.param("IUPAC")
        os.makedirs(dirname(to_file), exist_ok=True)
        cmd = r'''{_cat} {_file} | sed -r '/^[^>]/ {{ s/[{_IUPAC}]+/N/g; s/{_iupac}/n/g }}' > {_out}'''.format(
            _cat = self.is_gz(from_file) and "zcat" or "cat",
            _file = from_file,
            _IUPAC = IUPAC.upper(),
            _iupac = IUPAC.lower(),
            _out = to_file
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)


    # UTILS
    def is_gz(self, filename):
      return filename.endswith(".gz")


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

