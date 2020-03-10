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
  # out
  parser.add_argument("--meta_out", metavar="data/metadata/species", required = False,
                      type=argparse.FileType('w',  encoding='UTF-8'), default=sys.stdout,
                      help="stats output [STDOUT]" )
  parser.add_argument("--data_out_dir", metavar="data/metadata", required = True, type=str,
                      help="dir to store files into" )
  parser.add_argument("--genome_conf", metavar="genome.json", required = True, type=str,
                      help="genome json file output" )
  parser.add_argument("--manifest_out", metavar="manifest.json", required = True, type=str,
                      help="manifest file output" )
  # meta_defaults
  parser.add_argument("--species_division", metavar="EnsemblMetazoa", required = False,
                      type=str, default = "EnsemblMetazoa", help="species.division default" )
  parser.add_argument("--genebuild_method", metavar="import", required = False,
                      type=str, default = "import", help="genebuild.method default" )
  parser.add_argument("--genebuild_level", metavar="toplevel", required = False,
                      type=str, default = "toplevel", help="genebuild.level default" )
  #
  args = parser.parse_args()
  return args


## META CONF ##
class MetaConf:
  def __init__(self, config = None):
    self.tech_data = defaultdict(list)
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

  def dump(self, out):
    for k, vals in sorted(self.tech_data.items(), key=lambda x: x[0]):
      for v in sorted(vals):
        print("\t".join(["#CONF", str(k), str(v)]), file=out)
    for k, vals in sorted(self.data.items(), key=lambda x: x[0]):
      for v in vals:
        print("\t".join([str(k), str(v)]), file=out)

  def get(self, key, idx = 0, tech = False):
    d = self.data
    if tech:
      d = self.tech_data
    if key not in d or len(d[key]) < 1:
      return None
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
      if "db_xref" in qualifiers:
        taxon_id_pre = list(filter(lambda x: x.startswith("taxon:"), qualifiers["db_xref"]))[0]
        if taxon_id_pre:
          taxon_id = int(taxon_id_pre.split(":")[1])
          self.update("TAXON_ID", taxon_id, tech = True)

  def update_from_dict(self, d, k):
    if d is None:
      return
    if k not in d:
      return
    if not str(d[k]).strip():
      return
    self.update(k, d[k])

  def update_derived_data(self, defaults = None):
    # assembly metadata
    new_name = self.get("assembly.accession")
    if new_name:
      new_name = new_name.strip().replace("_","").replace(".","v")
      self.update("assembly.name", new_name)
    aname = self.get("assembly.name")
    self.update("assembly.default", aname)
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
      self.update("species.alias", [])
    # genebuild metadata
    self.update_from_dict(defaults, "genebuild.method")
    self.update_from_dict(defaults, "genebuild.level")
    self.update("genebuild.version", new_name + ".0")
    today = datetime.datetime.today()
    self.update("genebuild.start_date",
                "%s-%0d-%s" % (today.year, today.month, self.get("species.division")))
   


## GENOME CONF ##
class GenomeConf:
  def __init__(self, meta = None):
    pass
  def dump(self, out = None):
    pass


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
    shutil.copyfile(name, outfile)
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
    "species.division" : args.species_division,
    "genebuild.method" : args.genebuild_method,
    "genebuild.level" : args.genebuild_level,
  })
  meta.dump(args.meta_out)

  genome = GenomeConf(meta)
  genome.dump(args.genome_conf)

  manifest = Manifest({
    "fasta_dna" : args.fasta_dna,
    "fasta_pep" : args.fasta_pep,
    # "genome" : args.genome_conf.name, # check overriding on self copy (compare abs paths/ inodes?)
  })
  manifest.dump(args.manifest_out, outdir=args.data_out_dir)


if __name__ == "__main__":
    main()

