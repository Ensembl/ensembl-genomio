import dataclasses as dc
import gzip
import json
import os
import re
import sys

from collections import defaultdict, OrderedDict
from dataclasses import dataclass, field
from os.path import dirname, join as pj
from typing import Dict, List, Optional

from Bio import SeqIO # type: ignore

@dataclass
class SeqRegionSyn:
  name: str
  source: str

@dataclass
class SeqRegion:
  name: str
  length: int
  coord_system_level: str = "contig"
  location: Optional[str] = None
  circular : bool = False
  codon_table : Optional[int] = None
  synonyms: List[SeqRegionSyn] = field(default_factory=list)
  _rank: Optional[int] = None

def no_nulls_dict(x):
  return {k:v for k, v in x if v is not None}

class SeqRegionConf:
  def __init__(self,
               fasta_file: Optional[str],
               asm_rep_file: Optional[str],
               seq_region_raw: Optional[str],
               seq_region_genbank: Optional[str],
               meta: Optional[dict],
               syns_src: str = "GenBank"
               ):
    self.seq_regions : Dict[str, SeqRegion] = defaultdict(lambda:SeqRegion("",0))
    self.ordered : List[SeqRegion] = field(default_factory=list) # try like this
    self.syns : Dict[str, str] = dict()
    self.syn_src_default = syns_src
    # process
    #  get actual seq_region names
    self.fill_length_from_fasta(fasta_file)
    #  load data from seq_region_raw
    self.fill_info_from_seq_region_raw(seq_region_raw, use_syns = False)
    #  get data from meta
    self.fill_info_from_tech(meta)
    #  fill synonyms from report file
    self.fill_info_from_asm_rep(asm_rep_file)
    #  load from seq_region_genebank, infer names based on the loaded synonyms
    self.fill_info_from_seq_region_raw(seq_region_genbank, add_missing = False, update_from_new = True, update_syns = False)

  def dump(self, json_out: str) -> None:
    if not json_out:
      return
    if not self.seq_regions:
      return
    # reorder
    os.makedirs(dirname(json_out), exist_ok=True)
    with open(json_out, 'wt') as jf:
      out_list = list(map(lambda x: dc.asdict(x, dict_factory = no_nulls_dict), self.seq_regions.values()))
      # reorder based on _rank value and nullify it
      json.dump(out_list, jf, indent = 2, sort_keys = True)

  def merge_syns(self,
                 old_syns: List[SeqRegionSyn],
                 new_syns: Optional[List[SeqRegionSyn]],
                 use_new_source: bool = False) -> List[SeqRegionSyn]:
    if not new_syns:
      return old_syns
    if not old_syns:
      return new_syns
    syns_dict = defaultdict(set)
    for s in old_syns:
      syns_dict[s.name].add(s.source)
    for n in new_syns:
      if n.name not in syns_dict or not use_new_source:
        syns_dict[n.name].add(n.source)
      else:
        syns_dict[n.name]=set([n.source])
    out = []
    for name, sources in syns_dict.items():
      for src in sources:
        out.append(SeqRegionSyn(name, src))
    return out

  def fill_length_from_fasta(self, fasta_file: Optional[str]) -> None:
    """load actual seq_region nmaes and length from the fasta file if provided"""
    if not fasta_file:
      return
    _open = fasta_file.endswith(".gz") and gzip.open or open
    with _open(fasta_file, 'rt') as fasta:
      fasta_parser = SeqIO.parse(fasta, "fasta")
      for rec in fasta_parser:
        self.seq_regions[rec.name] = dc.replace(self.seq_regions[rec.name],
                                                name = rec.name,
                                                length = len(rec))
    return

  def fill_info_from_tech(self, tech_data: Optional[dict]) -> None:
    """load data synomyms from meta tags like CONTIG_CHR_..."""
    if not tech_data:
      return
    tk = tech_data.keys()
    # name for contig tags
    contig_chr_keys = list(filter(lambda x: x.upper().startswith("CONTIG_CHR_"), tk))

    # mito related keys
    mt_keys = frozenset(filter(lambda x: x.upper().startswith("MT_"), tk))
    mt_codon_table = int(tech_data.get("MT_CODON_TABLE", 0)) or None
    mt_circular = tech_data.get("MT_CIRCULAR", "0").upper() in ["YES", "1", "TRUE"]

    # default cs name for listed contigs
    cs_tag = tech_data.get("ORDERED_CS_TAG", "chromosome")

    for rank, k in enumerate(contig_chr_keys):
      syn = k.upper().split("_",2)[2]
      contig, *additional_syns = tech_data.get(k, "").split()
      # syns
      old_syns = self.seq_regions[contig].synonyms
      new_syns = [ SeqRegionSyn(s, s == syn and "Ensembl_Metazoa" or self.syn_src_default)
                    for s in ([syn] + additional_syns) ]
      syns = self.merge_syns(old_syns, new_syns, use_new_source = True)
      # update seq region
      self.seq_regions[contig] = dc.replace(
        self.seq_regions[contig],
        _rank = rank,
        name = contig,
        coord_system_level = cs_tag,
        synonyms = syns
      )
      # mito
      if mt_keys and syn == "MT":
        self.seq_regions[contig] = dc.replace(
          self.seq_regions[contig],
          _rank = rank,
          location = "mitochondrial_chromosome",
          circular = mt_circular,
          codon_table = mt_codon_table
        )
    return

  def fill_info_from_asm_rep(self, asm_rep_file: Optional[str], unversion = True) -> None:
    """load synonyms from assembly report, infer actual contig names,  add synonyms to them.
       should be called after loading from fasta or seq_region_raw"""
    if not asm_rep_file:
      return
    # INSDC_accession, INSDC_submitted_name
    use_cols = {
      "Sequence-Name" : "INSDC_submitted_name",
      "GenBank-Accn" : "INSDC",
      "RefSeq-Accn" : "RefSeq",
      "UCSC-style-name" : "GenBank",
    }

    _open = asm_rep_file.endswith(".gz") and gzip.open or open
    with _open(asm_rep_file, 'rt') as asm_rep:
      header_line = None
      header_fixed = None
      # 
      for line in asm_rep:
        if not header_fixed:
          if line.startswith("#"):
            header_line = line
            continue
          elif header_line:
            header_fixed = { i : use_cols.get(n.strip())
                               for i, n in enumerate(header_line[1:].split())
                                 if n.strip() in use_cols }
          else:
            break
        # split data line somehow
        data_fields = line.split("\t")
        if len(data_fields) <= 1:
          data_fields = line.split()

        sources_syns_raw = { src: data_fields[i].strip() for i, src in header_fixed.items() if i < len(data_fields) }
        sources_syns = { src: nm for src, nm in sources_syns_raw.items() if nm and nm.lower() != "na" }
        # remove INSDC_submitted_name if it's equal to INSDC or RefSeq one
        major_syns = { src: nm for src, nm in sources_syns.items() if src in ["INSDC", "RefSeq"] }
        if "INSDC_submitted_name" in sources_syns:
          if sources_syns["INSDC_submitted_name"] in set(major_syns.values()):
            sources_syns.pop("INSDC_submitted_name")

        # attach
        for src, syn in sources_syns.items():
          if src not in ["INSDC_submitted_name", "INSDC", "RefSeq"]:
            continue
          contig = syn
          if contig not in self.seq_regions:
            if unversion:
              contig = self.unversion(contig)
              if contig not in self.seq_regions:
                continue
              self.syns.update({syn : contig})
            else:
              continue
          self.syns.update({s : contig for s in sources_syns.values()})
          # update_syns
          old_syns = self.seq_regions[contig].synonyms
          new_syns = [ SeqRegionSyn(name, src) for src, name in sources_syns.items() if (name and src) ]
          syns = self.merge_syns(old_syns, new_syns, use_new_source = True)
          self.seq_regions[contig] = dc.replace(
            self.seq_regions[contig],
            synonyms = syns
          )
    return

  def unversion(self, name: str) -> str:
    if "." not in name:
      return name
    return re.sub(r'([^\.])\.\d+$', r'\1', name)

  def fill_info_from_seq_region_raw(self,
                                    raw_json: Optional[str],
                                    add_missing: bool = True,
                                    use_syns: bool = True,
                                    update_from_new: bool = True,
                                    update_syns: bool = True) -> None:
    if not raw_json:
      return
    _open = raw_json.endswith(".gz") and gzip.open or open
    with _open(raw_json, 'rt') as sr_raw:
      print("adding data from seq_regions json:", raw_json, file=sys.stderr)
      data = json.load(sr_raw, object_pairs_hook=OrderedDict)
      if type(data) != list:
        data = [ data ]
      for sr in data:
        name = sr.get("name")
        if not name:
          continue
        contig = name
        if use_syns and contig not in self.seq_regions:
          contig = self.syns.get(contig, contig)
        if not add_missing and contig not in self.seq_regions:
          continue
        # updating othewise
        sr_dict = dc.asdict(self.seq_regions[contig], dict_factory = no_nulls_dict)
        updates = {}
        for tag in "length coord_system_level location circular codon_table".split():
          val = sr.get(tag)
          if val is None: continue
          if not update_from_new and sr_dict.get(tag) is not None: continue
          updates[tag] = val
        # synonyms
        raw_syns = sr.get("synonyms")
        if update_syns and raw_syns:
          old_syns = self.seq_regions[contig].synonyms
          raw_syns = [ SeqRegionSyn(s["name"], s.get("source", self.syn_src_default)) for s in raw_syns ]
          updates["synonyms"] = self.merge_syns(old_syns, raw_syns, use_new_source = True)
        # update
        self.seq_regions[contig] = dc.replace(
          self.seq_regions[contig],
          **updates
        )
    return

