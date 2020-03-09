import argparse

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
  parser.add_argument("--data_out_dir", metavar="data/metadata", required = False, type=str,
                      help="dir to store files into" )
  parser.add_argument("--genome_conf", metavar="genome.json", required = True,
                      type=argparse.FileType('w',  encoding='UTF-8'),
                      help="genome json file output" )
  parser.add_argument("--manifest_out", metavar="manifest.json", required = True,
                      type=argparse.FileType('w',  encoding='UTF-8'),
                      help="manifest file output" )
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
      out[tag.strip()].append(rest[0])

  def dump(self, out):
    for k, vals in self.tech_data.items():
      for v in vals:
        print("\t".join(["#CONF", str(k), str(v)]), file=out)
    for k, vals in self.data.items():
      for v in vals:
        print("\t".join([str(k), str(v)]), file=out)

  def merge_from_gbff(self, gbff):
    if not gbff:
      return

    data = {}
    print("adding data from  %s" % gbff, file = sys.stderr)
    _open = gbff.endswith(".gz") and gzip.open or open
    with _open(fname, 'rt') as gbff:
      gb_parser = SeqIO.parse(gbff, "genbank")
      record = next(gb_parser)
      qualifiers = record.features[0].qualifiers

      data["scientific_name"] = qualifiers["organism"][0]
      data["strain"] = qualifiers["strain"][0]
      # strain.replace(" ","")?
      # production name
      taxon_id_pre = list(filter(lambda x: x.startswith("taxon:"), qualifiers["db_xref"]))[0]
      data["taxonomy_id"] = int(taxon_id_pre.split(":")[1])


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
    pass

  def is_gz(self, name):
    return name.endswith(".gz")

  def gunzip(self, name, tag, outdir):
    if not outdir:
      print("no out_dir specified to uncompress to", file=sys.stderr)
      return None

    sfx_pos_pre = name.rfind(".")
    sfx = name[name.rfind(".", 0, sfx_pos_pre):-1]
    outfile = pj(outdir,tag+sfx)

    os.makedirs(outdir, exist_ok=True)
    sp.run(r'''gunzip -c %s > %s''' % (name, outfile))
    return outfile

  def copy(self, name, tag, outdir):
    if not outdir:
      print("no out_dir specified to copy to", file=sys.stderr)
      return None

    sfx = name[name.rfind("."):-1]
    outfile = pj(outdir,tag+sfx)

    os.makedirs(outdir, exist_ok=True)
    shutil.copyfile(name, outfile)
    return outfile

  def md5sum(self, name):
    if not name:
      return
    pre = sp.check_output(["md5sum", name])
    (md5sum, *rest) = pre.split()
    return md5sum

  def dump(self, outfile, outdir = None):
    if outdir:
      os.makedirs(outdir, exist_ok=True)
    out = {}
    for tag, name in self.files.items():
      out_file = name
      if self.is_gz(name):
        out_file = self.gunzip(name, tag, outdir)
      elif self.always_copy:
        out_file = self.copy(name, tag, )
      md5 = md5sum(out_file)
      if out_file and md5:
        out[tag] = { "file" : out_file , "md5sum" : md5 }
    if out:
      json.dump(out, outfile, indent = 2)


## MAIN ##
def main():
  args = get_args()

  meta = MetaConf(args.raw_meta_conf)
  meta.merge_from_gbff(args.gbff_file)
  meta.dump(args.meta_out)

  genome = GenomeConf(meta)
  genome.dump(args.genome_conf)

  manifest = Manifest({
    "fasta_dna" : args.fasta_dna,
    "fasta_pep" : args.fasta_dna,
    "genome" : genome_out,
  })
  manifest.dump(args.manifest_out, outdir=args.data_out_dir)


if __name__ == "__main__":
    main()

