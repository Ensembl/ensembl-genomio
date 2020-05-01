import argparse
import datetime
import gzip
import json
import os
import re
import shutil
import subprocess as sp
import sys

from collections import defaultdict
from os.path import dirname, join as pj

from Bio import SeqIO


def get_args():
  parser = argparse.ArgumentParser()
  # in files
  parser.add_argument("--raw_meta_conf", metavar="meta/species", required=False,
                      type=argparse.FileType('rt', encoding='UTF-8'),
                      help="raw meta file with default settings")
  parser.add_argument("--gbff_file", metavar="species.gbff", required = False, type=str,
                      help="gbff file to get species name, accession, taxid from")
  parser.add_argument("--fasta_dna", metavar="species_dna.fna", required = True, type=str,
                      help="fasta file with dna sequences")
  parser.add_argument("--fasta_pep", metavar="species_pep.fa", required = False, type=str,
                      help="fasta file with protein sequences")
  parser.add_argument("--gff_file", metavar="species.gff3", required = False, type=str,
                      help="gff file with  gene models")
  parser.add_argument("--fann_file", metavar="functional_annotation.json", required = False, type=str,
                      help="json file with functional annotation")
  parser.add_argument("--seq_region_raw", metavar="seq_region_raw.json", required = False, type=str,
                      help="seq_region raw json file to patch")
  parser.add_argument("--asm_rep_file", metavar="species_assembly_report.txt", required = False, type=str,
                      help="GenBank assembly report to get seq_region syns from")
  # out
  parser.add_argument("--meta_out", metavar="data/metadata/species", required = False,
                      type=argparse.FileType('w',  encoding='UTF-8'), default=sys.stdout,
                      help="stats output [STDOUT]" )
  parser.add_argument("--data_out_dir", metavar="data/metadata", required = True, type=str,
                      help="dir to store files into" )
  parser.add_argument("--genome_conf", metavar="genome.json", required = True, type=str,
                      help="genome json file output" )
  parser.add_argument("--seq_region_conf", metavar="seq_region.json", required = False, type=str,
                      help="seq_region json file output" )
  parser.add_argument("--manifest_out", metavar="manifest.json", required = True, type=str,
                      help="manifest file output" )
  # meta_defaults
  parser.add_argument("--assembly_version", metavar="1", required = False,
                      type=int, default = 1, help="assembly.version default")
  parser.add_argument("--species_division", metavar="EnsemblMetazoa", required = False,
                      type=str, default = "EnsemblMetazoa", help="species.division default")
  parser.add_argument("--genebuild_method", metavar="import", required = False,
                      type=str, default = "import", help="genebuild.method default")
  parser.add_argument("--genebuild_level", metavar="toplevel", required = False,
                      type=str, default = "toplevel", help="genebuild.level default")
  parser.add_argument("--syns_src", metavar="GenBank", required = False,
                      type=str, default = "GenBank", help="syns source default")
  #
  args = parser.parse_args()
  return args


## META CONF ##
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
    _acc = self.get("assembly.accession").split(".")[0].replace("_", "")
    _strain = self.get("species.strain")
    _prod_name = ("%s_%s" % (_sci_name, _acc)).lower().replace(" ","_")
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
      #"annotation.provider_name",
      #"annotation.provider_url",
      #"assembly.provider_name",
      #"assembly.provider_url",
      "assembly.accession",
      "assembly.name",
      "genebuild.method",
      "genebuild.start_date",
      "genebuild.version",
      "provider.name",
      "provider.url",
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
                          syns_src = "GenBank"):
    if not json_out:
      return

    tk = self.tech_data.keys()
    chr_k = list(filter(lambda x: x.upper().startswith("CONTIG_CHR_"), tk))
    mt_k = frozenset(filter(lambda x: x.upper().startswith("MT_"), tk))

    ctg_len = dict()
    if fasta_file:
      _open = fasta_file.endswith(".gz") and gzip.open or open
      with _open(fasta_file, 'rt') as fasta:
        fasta_parser = SeqIO.parse(fasta, "fasta")
        for rec in fasta_parser:
          ctg_len[rec.name] = len(rec)

    asm_rep_syns = dict()
    if asm_rep_file:
      # INSDC_accession, INSDC_submitted_name 
      use_cols = {
        "Sequence-Name" : "INSDC_submitted_name",
        "GenBank-Accn" : "INSDC",
        "RefSeq-Accn" : "RefSeq",
        "Assigned-Molecule" : "GenBank",
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
          _data = line.split()
          _syns = { src: _data[i].strip() for i, src in header_fixed.items() if i < len(_data) }
          _out = { src: nm for src, nm in _syns.items() if nm and nm.lower() != "na" }
          for src in [ "INSDC", "RefSeq" ]: 
            if src not in _out:
              continue
            asm_rep_syns[_out[src]] = _out

    out = []
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
      out.append({
        "name" : ctg_id,
        "synonyms" : syns_out,
        "coord_system_level" : cs_tag,
      })
      if ctg_id in ctg_len:
        out[-1]["length"] = ctg_len[ctg_id]
      if mt_k and syn == "MT":
        out[-1]["location"] = "mitochondrial_chromosome"
        if "MT_CODON_TABLE" in mt_k:
          out[-1]["codon_table"] = int(self.get("MT_CODON_TABLE", tech=True))
        if "MT_CIRCULAR" in mt_k:
          mtc = self.get("MT_CIRCULAR", tech=True).strip().upper()
          out[-1]["circular"] = ( mtc == "YES" or mtc == "1" )

    used_ctg_names = frozenset([s["name"] for s in out])
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
      out.append(sr)

    # merge with seq_region_raw 
    if seq_region_raw:
      # load and merge into out
      pass

    if out:
      os.makedirs(dirname(json_out), exist_ok=True)
      with open(json_out, 'wt') as jf:
        json.dump(out, jf, indent = 2)


## MANIFEST CONF ##
class Manifest:
  def __init__(self, mapping, ungzip = True, always_copy = True):
    self.files = mapping
    self.ungzip = ungzip
    self.always_copy = always_copy

  def is_gz(self, name):
    return name.endswith(".gz")

  def gunzip(self, name, tag, outdir):
    if not outdir:
      print("no out_dir specified to uncompress to", file=sys.stderr)
      return None

    nogzname = name.replace(".gz", "")
    sfx = nogzname[nogzname.rfind("."):]
    outfile = pj(outdir,tag+sfx)

    os.makedirs(outdir, exist_ok=True)
    sp.run(r'''gunzip -c %s > %s''' % (name, outfile), shell = True)
    return outfile

  def copy(self, name, tag, outdir):
    if not outdir:
      print("no out_dir specified to copy to", file=sys.stderr)
      return None

    sfx = name[name.rfind("."):]
    outfile = pj(outdir,tag+sfx)

    os.makedirs(outdir, exist_ok=True)
    try:
      shutil.copyfile(name, outfile)
    except(shutil.SameFileError):
      pass
    return outfile

  def md5sum(self, name):
    if not name:
      return
    pre = sp.check_output("md5sum %s" % name, shell=True)
    (md5sum, *rest) = pre.split()
    return md5sum

  def dump(self, json_out, outdir = None):
    if outdir:
      os.makedirs(outdir, exist_ok=True)
    out = {}
    for tag, name in self.files.items():
      if not name:
        continue
      outfile = name
      if self.is_gz(name) and self.ungzip:
        print("uncompressing %s for %s" %(name, tag), file=sys.stderr)
        outfile = self.gunzip(name, tag, outdir)
      elif self.always_copy:
        print("copying %s for %s" %(name, tag), file=sys.stderr)
        outfile = self.copy(name, tag, outdir)
      print("calculating md5 for %s (%s)" %(outfile, tag), file=sys.stderr)
      md5 = self.md5sum(outfile)
      if outfile and md5:
        out[tag] = { "file" : outfile , "md5sum" : md5.decode() }
    if out:
      with open(json_out, 'wt') as jf:
        json.dump(out, jf, indent = 2)


## MAIN ##
def main():
  args = get_args()

  meta = MetaConf(args.raw_meta_conf)
  meta.merge_from_gbff(args.gbff_file)
  meta.update_derived_data({
    "assembly.version" : args.assembly_version,
    "genebuild.method" : args.genebuild_method,
    "genebuild.level" : args.genebuild_level,
    "species.division" : args.species_division,
  }, update_annotation_related = (args.gff_file is not None))
  meta.dump(args.meta_out)

  meta.dump_genome_conf(args.genome_conf)
  meta.dump_seq_region_conf(args.seq_region_conf,
                            fasta_file = args.fasta_dna,
                            asm_rep_file = args.asm_rep_file,
                            seq_region_raw = args.seq_region_raw,
                            syns_src = args.syns_src)

  manifest = Manifest({
    "fasta_dna" : args.fasta_dna,
    "fasta_pep" : args.fasta_pep,
    "genome" : args.genome_conf,
    "seq_region" : args.seq_region_conf,
    "functional_annotation" : args.fann_file,
    "gff3" : args.gff_file,
  })
  manifest.dump(args.manifest_out, outdir=args.data_out_dir)


if __name__ == "__main__":
    main()

