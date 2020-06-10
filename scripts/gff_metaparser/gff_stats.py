import argparse
import io
import json
import sys

from BCBio.GFF import GFFExaminer
from BCBio import GFF
from collections import defaultdict

# locals
from gffstruct.prefixtree import PfxTr
from gffstruct.validstruct import ValidStructures
from gffstruct.walkcontext import WalkContext


def get_args():
  parser = argparse.ArgumentParser()
  # configs
  parser.add_argument("--conf", metavar="structures.conf", required=False,
                      type=argparse.FileType('rt', encoding='UTF-8'),
                      help="valid structures config file")
  # various options
  parser.add_argument("--stats_only", action="store_true", required=False, help="produce only stats output")
  parser.add_argument("--fail_unknown", action="store_true", required=False, help="fail if unknown structure met")
  parser.add_argument("--processed_qualifiers",
                      metavar = "source,name,parent,dbxref,phase,product,protein_id,biotype",
                      required=False, type = str, default = "",
                      help="comma separated list of qualifiers to use [keeping all by default]")
  # ? separate sources, value to include #ALL_SOURCES with #
  parser.add_argument("--no_longest_pfx", action="store_true", required=False, help="no lognest pfrefixes for IDs" )
  parser.add_argument("--dump_used_options", action="store_true", required=False, help="dump used (not None) options")
  parser.add_argument("--rule_options", metavar="load_pseudogene_with_CDS,option2,...",
                     required=False, action="append",
                     help="options to control conditional rules")
  # output
  parser.add_argument("--stats_out", metavar="stats.out", required = False,
                      type=argparse.FileType('w',  encoding='UTF-8'), default=sys.stdout,
                      help="stats output [STDOUT]" )
  parser.add_argument("--detailed_report", metavar="detailed.log", required = False,
                      type=str, default=None,
                      help="file to report every failed and modified model" )
  parser.add_argument("--gff_out", metavar="out.gff3", required = False,
                      type=argparse.FileType('w', encoding='UTF-8'), default=sys.stdout,
                      help="resulting gff output [STDOUT]" )
  # input
  parser.add_argument("gff_in", metavar="in.gff3",
                      type=argparse.FileType('rt', encoding='UTF-8'),
                      help="input gff file (use '-' to read from STDIN)" )
  args = parser.parse_args()
  if args.dump_used_options:
    print("used options:", file=sys.stderr)
    for o, v in args.__dict__.items():
      if v: print("  %s: %s" %(o, v), file=sys.stderr)
  return args


def main():
  args = get_args()

  rule_options = args.rule_options
  if rule_options:
    rule_options = list(filter(None, map(lambda s: s.strip(), ",".join(rule_options).split(","))))
  known_structures = ValidStructures(args.conf, rule_options = rule_options)

  print("examining  %s" % (args.gff_in.name), file = sys.stderr)

  useful_qls = frozenset(filter(lambda x: x != "", args.processed_qualifiers.split(",")))
  print("using quals %s" % (args.processed_qualifiers), file = sys.stderr)

  stats = defaultdict(int)
  stats_sources = defaultdict(set)
  stats_prefixes = defaultdict(lambda: defaultdict(PfxTr)) # should be paths:item based

  def down_features(feature, prefix = None, prev_ids = []):
    #print(feature)
    # use https://docs.python.org/3/library/operator.html attrgetter methodcaller ??
    #print(feature.type, feature.id, feature.location)
    quals = feature.qualifiers
    out = { k : (len(v) > 0 and v[0] or None)
              for k,v in quals.items()
                if (not useful_qls or k.lower() in useful_qls) }
    #print(out)

    if prefix is None:
      prefix = feature.type
    else:
      prefix += "/" + feature.type
    new_prev_ids = prev_ids + [feature.id]

    if not feature.sub_features:
      known_structures.process(WalkContext(prefix), ignore_unseen = args.stats_only or args.conf is None)

      # analyze
      stats[prefix] += 1
      stats_sources[prefix].add(out["source"])
      for (ftype, fid) in zip (prefix.split("/"), new_prev_ids):
        stats_prefixes[prefix][ftype].add(fid)

    for sf in feature.sub_features:
      down_features(sf, prefix, new_prev_ids)


  gff =  GFF.parse(args.gff_in)
  for contig in gff:
    for cnt, feature in enumerate(contig.features):
      down_features(feature)
      # if cnt > 20: break
    # break


# create structures, provide link to top, and root ???
# set to None when done?

  for k in sorted(stats.keys(), key = lambda x: -stats[x]):
    pfx_info = {}
    if k in stats_prefixes:
      for part in stats_prefixes[k]:
        pfx_info[part] = stats_prefixes[k][part].get_max()
    print("%s\t%d\t%s\t%s" % (
      k, stats[k],
      ",".join(stats_sources[k]),
      json.dumps(pfx_info)
    ), file = args.stats_out)

    pass


if __name__ == "__main__":
    # execute only if run as a script
    main()


# zcat data/pfal/Pfalciparum.gff.gz | python new-genome-loader/scripts/gff_metaparser/gff_stats.py  --processed_qualifiers source,name,parent,dbxref,phase,product,protein_id,biotype --conf new-genome-loader/scripts/gff_metaparser/conf/valid_structures.conf -

# zcat data/pfal/Pfalciparum.gff.gz | python new-genome-loader/scripts/gff_metaparser/gff_stats.py  --stats-only --processed_qualifiers source,name,parent,dbxref,phase,product,protein_id,biotype -



