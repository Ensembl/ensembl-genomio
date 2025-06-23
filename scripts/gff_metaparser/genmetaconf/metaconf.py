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


## META CONF ##

import datetime
import gzip
import json
import re
import os
import sys

from collections import defaultdict
from os.path import abspath, dirname, join as pj
from typing import Dict, List, Optional

from Bio import SeqIO  # type: ignore

from .seqregionconf import SeqRegionConf


class MetaConf:
    def __init__(self, config=None, add_generated_species_aliases=False):
        self.tech_data = defaultdict(list)
        self._order = dict()
        self.data = defaultdict(list)

        self.add_generated_species_aliases = add_generated_species_aliases
        self.load_from_tsv(config)

    def load_from_tsv(self, tsv):
        if not tsv:
            return
        for raw in tsv:
            out = self.data
            is_tech = re.search(r"^\s*#\s*CONF\s+", raw)
            if is_tech:
                out = self.tech_data
                raw = raw[is_tech.span()[1] :]
            raw = raw.split("#")[0]
            if re.match(r"^\s*$", raw):
                continue
            tag, *rest = raw.split(maxsplit=1)
            if rest:
                out[tag.strip()].append(rest[0].rstrip())
                self._order[tag.strip()] = len(self._order)  # use the last rank for multi keys
        # adding raw file path
        meta_file_raw = abspath(tsv.name)
        self.tech_data["META_FILE_RAW"].append(meta_file_raw.rstrip())
        self._order["META_FILE_RAW"] = len(self._order)  # use the last rank for multi keys

    def dump(self, out):
        for k, vals in sorted(self.tech_data.items(), key=lambda x: x[0]):
            for v in sorted(vals):
                print("\t".join(["#CONF", str(k), str(v)]), file=out)
        for k, vals in sorted(self.data.items(), key=lambda x: x[0]):
            for v in vals:
                print("\t".join([str(k), str(v)]), file=out)

    def get(self, key, idx=0, tech=False, default=None):
        d = self.data
        if tech:
            d = self.tech_data
        if key not in d or len(d[key]) < 1:
            return default
        if idx is None:
            return d[key]
        if idx >= len(d[key]):
            return None
        return d[key][idx]

    def update(self, key, val, tech=False):
        out = self.data
        if tech:
            out = self.tech_data
        if val is None or not val:
            return
        if key in out and out[key]:
            return
        if isinstance(val, list):
            out[key] += val
        else:
            out[key].append(val)

    def normalise_asm_name(self, asm_name):
        if not asm_name:
            return asm_name
        asm_name = asm_name.strip()
        asm_name = re.sub(r"[^a-z0-9A-Z_\.]+", "_", asm_name)
        asm_name = re.sub(r"_+", "_", asm_name)
        asm_name = re.sub(r"^_+", "", asm_name)
        asm_name = re.sub(r"_+$", "", asm_name)
        return asm_name

    def merge_from_gbff(self, gbff_file):
        if not gbff_file:
            return

        print(f"adding data from {gbff_file}", file=sys.stderr)
        _open = gbff_file.endswith(".gz") and gzip.open or open
        with _open(gbff_file, "rt") as gbff:
            gb_parser = SeqIO.parse(gbff, "genbank")
            record = next(gb_parser)
            qualifiers = record.features[0].qualifiers
            if "organism" in qualifiers:
                sci_name = qualifiers["organism"][0]
                self.update("species.scientific_name", sci_name)
                self.update("organism.scientific_name", sci_name)
            if "strain" in qualifiers:
                strain = qualifiers["strain"][0]
                self.update("species.strain", strain)
                self.update("strain.type", "strain")
            elif "isolate" in qualifiers:
                strain = qualifiers["isolate"][0]
                self.update("species.strain", strain)
                self.update("strain.type", "isolate")
            elif "cultivar" in qualifiers:
                strain = qualifiers["cultivar"][0]
                self.update("species.strain", strain)
                self.update("strain.type", "cultivar")
            if "db_xref" in qualifiers:
                taxon_id_pre = list(filter(lambda x: x.startswith("taxon:"), qualifiers["db_xref"]))[0]
                if taxon_id_pre:
                    taxon_id = int(taxon_id_pre.split(":")[1])
                    self.update("TAXON_ID", taxon_id, tech=True)
            annotations = record.annotations
            if "structured_comment" in annotations:
                str_cmt = annotations["structured_comment"]
                if "Genome-Assembly-Data" in str_cmt:
                    gad = str_cmt["Genome-Assembly-Data"]
                    ankey = list(
                        filter(lambda x: "assembly" in x.lower() and "name" in x.lower(), gad.keys())
                    )
                    if ankey:
                        asm_name = gad[ankey[0]]
                        asm_name = self.normalise_asm_name(asm_name)
                        self.update("assembly.name", asm_name)

    def report_meta_value(self, line: str, pat: str) -> Optional[str]:
        if not line:
            return None
        if not pat:
            return None
        if re.match(pat, line):
            (_tag, value, *_rest) = line.split(sep=":", maxsplit=1)
            value = self.normalise_asm_name(value)
            return value
        return None

    def update_from_report_meta_value(self, line: str, pat: str, meta_key: str):
        if not line:
            return None
        if not pat:
            return None
        if not meta_key:
            return None
        value = self.report_meta_value(line, pat)
        if value:
            self.update(meta_key, value)

    def merge_from_asm_rep(self, asm_rep_file):
        if not asm_rep_file:
            return

        print(f"adding data from {asm_rep_file}", file=sys.stderr)
        _open = asm_rep_file.endswith(".gz") and gzip.open or open
        with _open(asm_rep_file, "rt") as asm_rep:
            for line in asm_rep:
                # Assembly name:  cgigas_uk_roslin_v1
                self.update_from_report_meta_value(line, r"#\s+Assembly name:", "assembly.name")
                # BioSample:      SAMEA110187692
                self.update_from_report_meta_value(line, r"#\s+BioSample:", "organism.biosample_id")
                # Assembly level: Chromosome
                self.update_from_report_meta_value(line, r"#\s+Assembly level:", "assembly.level")
                # GenBank assembly accession: GCA_947086385.1
                self.update_from_report_meta_value(
                    line, r"#\s+GenBank assembly accession:", "assembly.accession_insdc"
                )
                # RefSeq assembly accession: GCF_947086385.1
                self.update_from_report_meta_value(
                    line, r"#\s+RefSeq assembly accession:", "assembly.accession_refseq"
                )
                # RefSeq assembly and GenBank assemblies identical: yes

    def update_from_dict(self, d, k, tech=False):
        if d is None:
            return
        if k not in d:
            return
        if not str(d[k]).strip():
            return
        self.update(k, d[k], tech)

    def update_derived_data(self, defaults=None, update_annotation_related=False):
        # assembly metadata
        asm_acc = self.get("assembly.accession")
        new_name = str(asm_acc)
        if new_name:
            new_name = new_name.strip().replace("_", "").replace(".", "v")
            self.update("assembly.name", new_name)
        aname = self.get("assembly.name")
        self.update("assembly.default", aname)
        self.update_from_dict(defaults, "assembly.version", tech=True)
        # get annotation source
        ann_source = self.get("species.annotation_source", default="").strip()
        ann_source = self.normalise_asm_name(ann_source)
        # picking assembly.alt_accession
        asm_acc_insdc = self.get("assembly.accession_insdc")
        asm_acc_refseq = self.get("assembly.accession_refseq")
        if not self.get("assembly.alt_accession"):
            if asm_acc_refseq and asm_acc_insdc:
                # only if species.annotaion_source is ~ "RefSeq"
                if "refseq" in ann_source.lower():
                    if asm_acc.startswith("GCF_"):
                        self.update("assembly.alt_accession", asm_acc_insdc)
                    else:
                        self.update("assembly.alt_accession", asm_acc_refseq)
        asm_acc_main = asm_acc_insdc or asm_acc
        # species metadata
        _acc = str(asm_acc).replace("_", "").replace(".", "v")
        _sci_name = self.get("species.scientific_name", default="")
        _prod_name_pfx = _sci_name.strip().lower()
        _prod_name_pfx = "_".join(re.sub(r"[^a-z0-9A-Z]+", "_", _prod_name_pfx).split("_")[:2])
        _prod_name = ("%s_%s" % (_prod_name_pfx, _acc)).lower().replace(" ", "_")
        _acc_main = str(asm_acc_main).replace("_", "").replace(".", "v")
        _prod_name_main = ("%s_%s" % (_prod_name_pfx, _acc_main)).lower().replace(" ", "_")
        #
        _strain = self.get("species.strain")
        self.update_from_dict(defaults, "species.division")
        # possibly add annotation source suffix
        _ann_source_sfx = self.get("ANNOTATION_SOURCE_SFX", tech=True, default="").strip()
        _ann_source_sfx = self.normalise_asm_name(_ann_source_sfx).replace("_", "").lower()[:2]
        if _ann_source_sfx:
            self.update("ANNOTATION_SOURCE_SFX", _ann_source_sfx, tech=True)
            _prod_name += _ann_source_sfx
            _prod_name_main += _ann_source_sfx
        self.update("species.production_name", _prod_name)
        self.update("species.production_name_main", _prod_name_main)
        _comm_name = self.get("species.common_name")
        _display_name = _sci_name
        if _strain or _comm_name:
            _strain_comm_part = ", ".join(map(str, filter(None, [_comm_name, _strain])))
            _display_name += f" ({_strain_comm_part})"
        _display_name_main = _display_name
        _display_name += " - " + asm_acc
        # possibly add annotation source tag
        _ann_source = self.get("species.annotation_source", default="").strip()
        _ann_source = self.normalise_asm_name(_ann_source)
        if _ann_source:
            self.update("species.annotation_source", _ann_source)
            self.update("genebuild.annotation_source", _ann_source)
            _display_name = f"{_display_name} [{_ann_source} annotation]"
        self.update("species.display_name", _display_name)
        self.update("species.display_name_main", _display_name_main)
        # back to using "Binomial_name_GCA_000001.1rs" names, only for GenBank ('GCA') accessions
        # get a s/_gc([af])(\d+)v(\d+)/_GC\U\1_\2.\3/i equivalent
        gc_map = lambda m: f"_GC{m.group(1).upper()}_{m.group(2)}.{m.group(3)}"
        # url
        _species_url = _prod_name.capitalize()
        _species_url = re.sub(r"_gc([af])(\d+)v(\d+)", gc_map, _species_url, flags=re.I)
        self.update("species.url", _species_url)
        # url for the main
        _species_url_main = _prod_name_main.capitalize()
        _species_url_main = re.sub(r"_gc([af])(\d+)v(\d+)", gc_map, _species_url_main, flags=re.I)
        self.update("species.url_main", _species_url_main)
        # syns
        syns = []
        sci_name_words = list(filter(None, _sci_name.split()))
        if sci_name_words:
            syns.append(sci_name_words[0][0] + ". " + sci_name_words[1])
            syns.append(sci_name_words[0][0] + "." + sci_name_words[1][:5])
            syns.append((sci_name_words[0][0] + sci_name_words[1][:5]).lower())
        if syns and self.add_generated_species_aliases:
            self.update("species.alias", syns)
        # organism metadata
        # adding duplicates to deal with RapidRelease/MVP requirements
        taxon_id = int(self.get("TAXON_ID", tech=True))
        self.update("organism.taxonomy_id", taxon_id)
        self.update("organism.species_taxonomy_id", taxon_id)
        self.update("organism.strain", _strain)
        self.update("organism.strain_type", self.get("strain.type"))
        self.update("organism.production_name", _prod_name)
        self.update("organism.common_name", _comm_name)
        # genebuild metadata
        if update_annotation_related:
            self.update_from_dict(defaults, "genebuild.method")
            self.update_from_dict(defaults, "genebuild.method_display")
            self.update_from_dict(defaults, "genebuild.level")
            self.update("genebuild.version", aname.replace("_", "").replace(".", "v") + ".0")
            today = datetime.datetime.today()
            self.update(
                "genebuild.start_date", "%s-%02d-%s" % (today.year, today.month, self.get("species.division"))
            )
            self.update("genebuild.initial_release_date", "%s-%02d" % (today.year, today.month))
            self.update("genebuild.last_geneset_update", "%s-%02d" % (today.year, today.month))
            # MVP/RR duplicates
            self.update("genebuild.provider_name", self.get("annotation.provider_name"))
            self.update("genebuild.provider_url", self.get("annotation.provider_url"))

    def dump_genome_conf(self, json_out):
        out = {}
        fields = [
            "annotation.provider_name",
            "annotation.provider_url",
            "assembly.provider_name",
            "assembly.provider_url",
            "assembly.accession",
            # not yet supported in genome.json schema
            # "assembly.accession_insdc",
            # "assembly.accession_refseq",
            # "assembly.alt_accession",
            # "assembly.level",
            "assembly.name",
            "genebuild.method",
            "genebuild.method_display",
            # "genebuild.provider_name",
            # "genebuild.provider_url",
            # "genebuild.annotation_source",
            "genebuild.start_date",
            "genebuild.version",
            # "organism.biosample_id",
            # "organism.common_name",
            # "organism.ensembl_name",
            # "organism.scientific_name",
            # "organism.strain",
            # "organism.strain_type",
            # "organism.species_taxonomy_id",
            # "organism.taxonomy_id",
            "*species.alias",
            "species.annotation_source",
            "species.display_name",
            "species.division",
            "species.production_name",
            "species.scientific_name",
            "species.strain",
            # "strain.type",
        ]
        for f in fields:
            if f.startswith("*"):
                self.split_add(out, f[1:], self.get(f[1:], idx=None))
            else:
                self.split_add(out, f, self.get(f))
        self.split_add(out, "assembly.version", self.get("assembly.version", tech=True))
        self.split_add(out, "species.taxonomy_id", int(self.get("TAXON_ID", tech=True)))

        # get chr aliases
        tk = self.tech_data.keys()
        chr_k = list(filter(lambda x: x.upper().startswith("CONTIG_CHR_"), tk))
        if chr_k:
            ctg_lst = [self.get(k, tech=True).split()[0] for k in sorted(chr_k, key=lambda x: self._order[x])]
            out["assembly"]["chromosome_display_order"] = ctg_lst

        if out:
            os.makedirs(dirname(json_out), exist_ok=True)
            with open(json_out, "wt") as jf:
                json.dump(out, jf, indent=2)

    def split_add(self, out, key, val):
        if val is None:
            return
        keys = key.split(".")
        if not keys:
            return
        pre = {keys[-1]: val}
        for k in keys[:-1]:
            if k not in out:
                out[k] = dict()
            out = out[k]
        out.update(pre)

    def dump_seq_region_conf(
        self,
        json_out,
        fasta_file=None,
        asm_rep_file=None,
        seq_region_raw=None,
        seq_region_genbank=None,
        seq_region_syns=None,
        syns_src="GenBank",
        default_genetic_code=1,
        default_circular=False,
    ):
        if not json_out:
            return
        sr_conf = SeqRegionConf(
            fasta_file=fasta_file,
            asm_rep_file=asm_rep_file,
            seq_region_raw=seq_region_raw,
            seq_region_genbank=seq_region_genbank,
            seq_region_syns=seq_region_syns,
            syns_src=syns_src,
            meta=self.tech_data,
            default_genetic_code=default_genetic_code,
            default_circular=default_circular,
        )
        sr_conf.dump(json_out)
