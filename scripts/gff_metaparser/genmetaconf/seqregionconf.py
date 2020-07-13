import dataclasses as dc
import gzip
import json
import os

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
    #self.ordered : List[SeqRegion]
    self.syns : Dict[str, str]  = field(default_factory=dict)
    # process
    #  get actual seq_region names
    self.fill_length_from_fasta(fasta_file)
    return
    #  load data from seq_region_raw
    self.fill_info_from_seq_region_raw(seq_region_raw)
    #  get data from meta
    self.fill_info_from_tech(meta)
    #  fill synonyms from report file
    self.fill_info_from_asm_rep(asm_rep_file)
    #  load from seq_region_genebank, infer names based on the loaded synonyms
    self.fill_info_from_seq_region_raw(seq_region_genbank, ignore_syns = True)

  def dump(self, json_out: str) -> None:
    if not json_out:
      return
    if not self.seq_regions:
      return
    # reorder
    os.makedirs(dirname(json_out), exist_ok=True)
    with open(json_out, 'wt') as jf:
      out_list = list(map(lambda x: dc.asdict(x, dict_factory = no_nulls_dict), self.seq_regions.values()))
      json.dump(out_list, jf, indent = 2, sort_keys = True)

  def fill_length_from_fasta(self, fasta_file: Optional[str]) -> None:
    """load actual seq_region nmaes and length from the fasta file if provided"""
    if fasta_file:
      _open = fasta_file.endswith(".gz") and gzip.open or open
      with _open(fasta_file, 'rt') as fasta:
        fasta_parser = SeqIO.parse(fasta, "fasta")
        for rec in fasta_parser:
          self.seq_regions[rec.name] = dc.replace(self.seq_regions[rec.name],
                                                  name = rec.name,
                                                  length = len(rec))
    return

  def fill_info_from_tech(self, tech_data: Optional[dict]) -> None:
    tk = self.tech_data.keys()
    chr_k = list(filter(lambda x: x.upper().startswith("CONTIG_CHR_"), tk))
    mt_k = frozenset(filter(lambda x: x.upper().startswith("MT_"), tk))

    out = OrderedDict()
    for k in chr_k:
      ctg_id, *syns = self.get(k, tech = True).split()
      syn = k.upper().split("_", 2)[2]
      syns.append(syn)
      syns = list(set(syns))
      syns_out = { nm : (nm == syn and "Ensembl_Metazoa" or syns_src) for nm in syns }
      if ctg_id in asm_rep_syns:
        for src, nm in asm_rep_syns[ctg_id].items():
          syns_out[nm] = src
      syns_out = [ { "name" : nm, "source" : src } for nm, src in syns_out.items() ]
      # merge syns_out with asm_rep_syns
      cs_tag = self.get("ORDERED_CS_TAG", tech = True, default = "chromosome")
      out[ctg_id] = {
        "name" : ctg_id,
        "synonyms" : syns_out,
        "coord_system_level" : cs_tag,
      }
      if ctg_id in ctg_len:
        out[ctg_id]["length"] = ctg_len[ctg_id]
      if mt_k and syn == "MT":
        out[ctg_id]["location"] = "mitochondrial_chromosome"
        if "MT_CODON_TABLE" in mt_k:
          out[ctg_id]["codon_table"] = int(self.get("MT_CODON_TABLE", tech=True))
        if "MT_CIRCULAR" in mt_k:
          mtc = self.get("MT_CIRCULAR", tech=True).strip().upper()
          out[ctg_id]["circular"] = ( mtc == "YES" or mtc == "1" )

    used_ctg_names = frozenset(out.keys())
    for ctg_id in ctg_len:
      if ctg_id in used_ctg_names:
        continue
      sr = {
        "name" : ctg_id,
        "length" : ctg_len[ctg_id],
        "coord_system_level" : "contig",
      }
      if ctg_id in asm_rep_syns:
        syns_out = [ { "name" : nm, "source" : src } for src, nm in asm_rep_syns[ctg_id].items() ]
        if syns_out:
          sr["synonyms"] = syns_out
      out[ctg_id] = sr


  def fill_info_from_asm_rep(self) -> None:
    asm_rep_syns = dict()
    if asm_rep_file:
      # INSDC_accession, INSDC_submitted_name
      use_cols = {
        "Sequence-Name" : "INSDC_submitted_name",
        "GenBank-Accn" : "INSDC",
        "RefSeq-Accn" : "RefSeq",
     }
      _open = asm_rep_file.endswith(".gz") and gzip.open or open
      with _open(asm_rep_file, 'rt') as asm_rep:
        header_line = None
        header_fixed = None
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
          _data = line.split("\t")
          if len(_data) <= 1:
            _data = line.split()
          _syns = { src: _data[i].strip() for i, src in header_fixed.items() if i < len(_data) }
          _out = { src: nm for src, nm in _syns.items() if nm and nm.lower() != "na" }
          if "INSDC_submitted_name" in _out and "INSDC" in _out and _out["INSDC_submitted_name"] == _out["INSDC"]:
            _out.pop("INSDC_submitted_name")
          if "INSDC_submitted_name" in _out and "RefSeq" in _out and _out["INSDC_submitted_name"] == _out["RefSeq"]:
            _out.pop("INSDC_submitted_name")
          for src in [ "INSDC", "RefSeq" ]:
            if src not in _out:
              continue
            asm_rep_syns[_out[src]] = _out
    return

  def fill_info_from_seq_region_raw(self, raw_json: str, ignore_syns: bool = False) -> None:
    # merge with seq_region_raw
    out_list = []
    merged_seq_regions = set()
    if seq_region_raw:
      _open = seq_region_raw.endswith(".gz") and gzip.open or open
      with _open(seq_region_raw, 'rt') as sr_raw:
        _data = json.load(sr_raw, object_pairs_hook=OrderedDict)
        if type(_data) != list:
          _data = [ _data ]
        for rawsr in _data:
          if "name" not in rawsr or not rawsr["name"]:
            continue
          name = rawsr["name"]
          # update synonym sources
          for syn in filter(lambda s: "source" not in s, rawsr.get("synonyms", [])):
            syn["source"] = syns_src
          # merge with out data
          merged_seq_regions.add(name)
          if name not in out:
            out_list.append(rawsr)
          else:
            rawsyns = rawsr.get("synonyms", [])
            out_data = out[name]
            for k in out_data:
              if k not in rawsr or (not rawsr[k] and out_data[k]):
                rawsr[k] = out_data[k]
            # more accurate merge of syns, not sure about "coord_system_level"
            if rawsyns:
              # put everything that left
              used_syns = frozenset([s["name"] for s in rawsyns])
              rawsyns += [ s for s in out_data.get("synonyms", []) if s["name"] not in used_syns ]
            # append
            out_list.append(rawsr)

    # add unused seq_regions from out
    for name, sr in out.items():
      if name not in merged_seq_regions:
        out_list.append(sr)


