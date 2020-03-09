import gzip
import json
import sys

from Bio import SeqIO

fname = sys.argv[1]
#opener = fname.endswith(".gz") and gzip.open or open
print("parsing %s" % fname, file = sys.stderr)

with gzip.open(fname, 'rt') as gbff:
  gb_parser = SeqIO.parse(gbff, "genbank")
  record = next(gb_parser)
  qualifiers = record.features[0].qualifiers
  data = {}
  data["scientific_name"] = qualifiers["organism"][0]
  data["strain"] = qualifiers["strain"][0]
  # strain.replace(" ","")?
  # production name
  taxon_id_pre = list(filter(lambda x: x.startswith("taxon:"), qualifiers["db_xref"]))[0]
  data["taxonomy_id"] = int(taxon_id_pre.split(":")[1])

  json.dump(data, sys.stdout, indent = 2)
