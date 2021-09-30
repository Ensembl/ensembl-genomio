## META CONF ##

import datetime
import gzip
import json
import re
import os
import sys

from collections import defaultdict
from os.path import abspath, dirname, join as pj

from Bio import SeqIO # type: ignore

from .seqregionconf import SeqRegionConf


class MetaConf:
  def __init__(self, config = None):
    self.tech_data = defaultdict(list)
    self._order = dict()
    self.data = defaultdict(list)

    self.load_from_tsv(config)

  def load_from_tsv(self, tsv):
    if not tsv:
      return
    for raw in tsv:
      out = self.data
      is_tech = re.search(r'^\s*#\s*CONF\s+', raw)
      if is_tech:
        out = self.tech_data
        raw = raw[is_tech.span()[1]:]
      raw = raw.split("#")[0]
      if re.match(r'^\s*$', raw):
        continue
      tag, *rest = raw.split(maxsplit = 1)
      if rest:
        out[tag.strip()].append(rest[0].rstrip())
        self._order[tag.strip()] = len(self._order) # use the last rank for multi keys
    # adding raw file path
    meta_file_raw = abspath(tsv.name)
    self.tech_data["META_FILE_RAW"].append(meta_file_raw.rstrip())
    self._order["META_FILE_RAW"] = len(self._order) # use the last rank for multi keys


  def dump(self, out):
    for k, vals in sorted(self.tech_data.items(), key=lambda x: x[0]):
      for v in sorted(vals):
        print("\t".join(["#CONF", str(k), str(v)]), file=out)
    for k, vals in sorted(self.data.items(), key=lambda x: x[0]):
      for v in vals:
        print("\t".join([str(k), str(v)]), file=out)

  def get(self, key, idx = 0, tech = False, default = None):
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

  def update(self, key, val, tech = False):
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

  def merge_from_gbff(self, gbff_file):
    if not gbff_file:
      return

    print("adding data from  %s" % gbff_file, file = sys.stderr)
    _open = gbff_file.endswith(".gz") and gzip.open or open
    with _open(gbff_file, 'rt') as gbff:
      gb_parser = SeqIO.parse(gbff, "genbank")
      record = next(gb_parser)
      qualifiers = record.features[0].qualifiers
      if "organism" in qualifiers:
        sci_name = qualifiers["organism"][0]
        self.update("species.scientific_name", sci_name)
      if "strain" in qualifiers:
        strain = qualifiers["strain"][0]
        self.update("species.strain", strain)
      elif "isolate" in qualifiers:
        strain = qualifiers["isolate"][0]
        self.update("species.strain", strain)
      if "db_xref" in qualifiers:
        taxon_id_pre = list(filter(lambda x: x.startswith("taxon:"), qualifiers["db_xref"]))[0]
        if taxon_id_pre:
          taxon_id = int(taxon_id_pre.split(":")[1])
          self.update("TAXON_ID", taxon_id, tech = True)
      annotations = record.annotations
      if "structured_comment" in annotations:
        str_cmt = annotations["structured_comment"]
        if "Genome-Assembly-Data" in str_cmt:
          gad = str_cmt["Genome-Assembly-Data"]
          ankey = list(filter(lambda x: "assembly" in x.lower() and "name" in x.lower(), gad.keys()))
          if ankey:
            self.update("assembly.name", gad[ankey[0]])

  def update_from_dict(self, d, k, tech = False):
    if d is None:
      return
    if k not in d:
      return
    if not str(d[k]).strip():
      return
    self.update(k, d[k], tech)

  def update_derived_data(self, defaults = None, update_annotation_related = False):
    # assembly metadata
    new_name = self.get("assembly.accession")
    if new_name:
      new_name = new_name.strip().replace("_","").replace(".","v")
      self.update("assembly.name", new_name)
    aname = self.get("assembly.name")
    self.update("assembly.default", aname)
    self.update_from_dict(defaults, "assembly.version", tech = True)
    # species metadata
    self.update_from_dict(defaults, "species.division")
    _sci_name = self.get("species.scientific_name")
    _acc = self.get("assembly.accession").replace("_","").replace(".","v")
    _strain = self.get("species.strain")
    _prod_name = _sci_name.strip().lower()
    _prod_name = "_".join(re.sub(r'[^a-z0-9A-Z]+', '_', _prod_name).split("_")[:2])
    _prod_name = ("%s_%s" % (_prod_name, _acc)).lower().replace(" ","_")
    self.update("species.production_name", _prod_name)
    _display_name = _sci_name
    if _strain:
      _display_name = ("%s (%s)" % (_sci_name, _strain))
    self.update("species.display_name", _display_name)
    self.update("species.url", _prod_name.capitalize())
    # syns
    syns = []
    w = list(filter(None, _sci_name.split()))
    syns.append(w[0][0] + ". " + w[1])
    syns.append(w[0][0] + "." + w[1][:3])
    syns.append((w[0][0] + w[1][:3]).lower())
    if syns:
      self.update("species.alias", syns)
    # genebuild metadata
    if update_annotation_related:
      self.update_from_dict(defaults, "genebuild.method")
      self.update_from_dict(defaults, "genebuild.level")
      self.update("genebuild.version", aname.replace("_", "").replace(".","v") + ".0")
      today = datetime.datetime.today()
      self.update("genebuild.start_date",
                  "%s-%02d-%s" % (today.year, today.month, self.get("species.division")))

  def dump_genome_conf(self, json_out):
    out = {}
    fields = [
      "annotation.provider_name",
      "annotation.provider_url",
      "assembly.provider_name",
      "assembly.provider_url",
      "assembly.accession",
      "assembly.name",
      "genebuild.method",
      "genebuild.start_date",
      "genebuild.version",
      "*species.alias",
      "species.display_name",
      "species.division",
      "species.production_name",
      "species.scientific_name",
      "species.strain",
    ]
    for f in fields:
      if f.startswith("*"):
        self.split_add(out, f[1:], self.get(f[1:], idx=None))
      else:
        self.split_add(out, f, self.get(f))
    self.split_add(out, "assembly.version", self.get("assembly.version", tech=True))
    self.split_add(out, "species.taxonomy_id", self.get("TAXON_ID", tech=True))

    # get chr aliases
    tk = self.tech_data.keys()
    chr_k = list(filter(lambda x: x.upper().startswith("CONTIG_CHR_"), tk))
    if chr_k:
       ctg_lst = [ self.get(k, tech = True).split()[0] for k in sorted(chr_k, key = lambda x: self._order[x]) ]
       out["assembly"]["chromosome_display_order"] = ctg_lst

    if out:
      os.makedirs(dirname(json_out), exist_ok=True)
      with open(json_out, 'wt') as jf:
        json.dump(out, jf, indent = 2)

  def split_add(self, out, key, val):
    if val is None:
      return
    keys = key.split(".")
    if not keys:
      return
    pre = {keys[-1]: val}
    for k in keys[:-1]:
      if k not in out:
        out[k] =  dict()
      out = out[k]
    out.update(pre)

  def dump_seq_region_conf(self, json_out,
                          fasta_file = None, asm_rep_file = None, seq_region_raw = None,
                          seq_region_genbank = None, seq_region_syns = None,
                          syns_src = "GenBank", default_genetic_code = 1):
    if not json_out:
      return
    sr_conf = SeqRegionConf(fasta_file = fasta_file,
                            asm_rep_file = asm_rep_file,
                            seq_region_raw = seq_region_raw,
                            seq_region_genbank = seq_region_genbank,
                            seq_region_syns = seq_region_syns,
                            syns_src = syns_src,
                            meta = self.tech_data,
                            default_genetic_code = default_genetic_code)
    sr_conf.dump(json_out)
