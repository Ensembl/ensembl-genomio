import io
import gzip
import json
import pprint
import sys

from BCBio.GFF import GFFExaminer
from BCBio import GFF
from collections import defaultdict


class PfxTr:
  def __init__(self):
    self.letters = dict()
    self.count = int(0)

  def add(self, s):
    if s is None or s == "":
      return
    root = self
    for c in [""] + list(s):
      if c not in root.letters:
        root.letters[c] = PfxTr()
      root = root.letters[c]
      root.count += 1

  def get_max(self):
     (pfx, max_cnt) = ("", 0)
     root = self
     while root.letters:
       c = max(root.letters.keys(), key = lambda x: root.letters[x].count)
       if root.letters[c].count < max_cnt:
         return(pfx, max_cnt)
       max_cnt = root.letters[c].count
       pfx += c
       root = root.letters[c] 
     return (pfx, max_cnt)


gff_file = sys.argv[1]
print("examining  %s" % (gff_file), file = sys.stderr)

#gff_types = {}
#with io.TextIOWrapper(gzip.open(sys.argv[1], 'r')) as in_handle:
#  examiner = GFFExaminer()
#  limits = examiner.available_limits(in_handle)
#  print("seen source and types: ", limits['gff_source_type'], file = sys.stderr)
#  print("seen types: ", limits['gff_type'], file = sys.stderr)
#  gff_types = { k[0].lower() : k[0] for k in limits['gff_type']}

useful_quals = frozenset("source Name Parent Dbxref phase product protein_id biotype".split())

stats = defaultdict(int)
stats_sources = defaultdict(set)
stats_prefixes = defaultdict(lambda: defaultdict(PfxTr)) # should be paths:item based

def down_features(feature, prefix = None, prev_ids = []):
   #print(feature)
   #print(feature.type, feature.id, feature.location)
   quals = feature.qualifiers
   out = { k:v[0] for k,v in quals.items() if k in useful_quals }
   #print(out)

   if prefix is None:
     prefix = feature.type
   else:
      prefix += "/" + feature.type
   new_prev_ids = prev_ids + [feature.id]

   if not feature.sub_features:
     # analyze
     stats[prefix] += 1
     stats_sources[prefix].add(out["source"])
     for (ftype, fid) in zip (prefix.split("/"), new_prev_ids):
        stats_prefixes[prefix][ftype].add(fid)

   for sf in feature.sub_features:
     down_features(sf, prefix, new_prev_ids)


with io.TextIOWrapper(gzip.open(sys.argv[1], 'r')) as in_handle:
  gff =  GFF.parse(in_handle)
  for contig in gff:
    for cnt, feature in enumerate(contig.features):
      down_features(feature)
      # if cnt > 3: break
    #break

for k in sorted(stats.keys(), key = lambda x: -stats[x]):
   pfx_info = {}
   if k in stats_prefixes:
     for part in stats_prefixes[k]:
       pfx_info[part] = stats_prefixes[k][part].get_max()
   print("%s\t%d\t%s\t%s" % (
     k, stats[k],
     ",".join(stats_sources[k]),
     json.dumps(pfx_info)
   ))


#  limits = examiner.available_limits(in_handle)
#  limits.keys()
#  #Out[9]: dict_keys(['gff_id', 'gff_source_type', 'gff_source', 'gff_type'])

#  limits['gff_source_type']
#  #{('Genbank', 'region'): 459,
# ('Genbank', 'gene'): 15577,
# ('Genbank', 'mRNA'): 15577,
# ('Genbank', 'exon'): 74199,
# ('Genbank', 'CDS'): 68827}

# list(limits['gff_source_type'].keys())[0][1]
# #Out[13]: 'region'

#limits['gff_type']
#Out[14]:
#{('region',): 459,
# ('gene',): 15577,
# ('mRNA',): 15577,
# ('exon',): 74199,
# ('CDS',): 68827}

#limits['gff_source']
#Out[15]: {('Genbank',): 174639}






