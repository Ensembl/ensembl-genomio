import functools
import re
import sys

from BCBio import GFF

# locals
from .seqfeatcomposer import SeqFeatComposer
from .utils import UpdatingLen
from .walkcontext import WalkContext


class GFF3Walker:
  _valid_structure_tags = frozenset("fullPath anyQual".split()) # /gene/mrna/exon vs exon/id, exon


  def __init__(self, parser, in_file, structure_tags="fullPath", global_ctx = None, norm_id = None):
    if not parser:
      raise Exception("not a valid parser: %s" % str(parser))
    if structure_tags not in self._valid_structure_tags:
      raise Exception("not a valid 'structure_tags' value: %s" % structure_tags)
    self._parser = parser
    self._struct_tags = structure_tags
    self._in_file = in_file
    self.global_context = global_ctx
    self._norm_id = norm_id
    self._supported_fields = dict()


  def norm_id(self, id, type=None):
    if self._norm_id is None:
      return id
    return self._norm_id(id, type = type)


  def supporting(self, obj, field):
    tp = type(obj)
    if tp not in self._supported_fields:
      self._supported_fields[tp] = frozenset(obj.__dir__())
    return field in self._supported_fields[tp]


  def walk(self, out_file = None, seq_len_dict = None):
    gff =  GFF.parse(self._in_file)
    sf_composer = SeqFeatComposer() #  ?? init composer once, get on __init__???
    for contig in gff:
      # get len from gff or fna file, try to infere if not available
      ctg_len = UpdatingLen(len(contig))
      ctg_len.update(seq_len_dict(contig.id), stop_on_success = True)
      # feature processing
      for cnt, topft in enumerate(contig.features):
        # iterate through all of the features
        context = WalkContext(global_context = self.global_context, ctg_len_inferer = ctg_len)
        context.update("_SEQID", contig.id)
        self.process_feature(topft, context)
      # dump gff
      sf_composer.gff_write(context, contig_len = len(ctg_len), out_file = out_file)


  def process_feature(self, feat, context):
    # store stats in the context
    loc = feat.location
    quals = feat.qualifiers
    phase = quals.get("phase")
    fulltag_parts = list(filter(None,[context.get("_FULLTAG"), feat.type]))
    fulltag = "/".join(fulltag_parts)
    depth = fulltag.count("/") + 1
    source = self.supporting(feat, "source") and feat.source or None
    is_leaf = not self.supporting(feat, "sub_features") or not feat.sub_features # fullPath processed and ID can be ommited only for leaves

    #norm
    if self._norm_id and self._norm_id.supporting(feat.type):
      feat.id = self.norm_id(feat.id, feat.type)
      if quals and quals.get("ID"):
        quals["ID"][0] = self.norm_id(quals["ID"][0], feat.type)

    if depth == 1:
      # set top feature
      context.top(feat)

    context.update(
      force_clean = True,
      _ID      = feat.id,
      _SRC     = source,
      _TYPE    = feat.type,
      _START   = loc.start,
      _END     = loc.end,
      _STRAND  = feat.strand,
      _PHASE   = phase,
      _QUALS   = quals,
      _FULLTAG = fulltag,
      _DEPTH   = depth,
      _ISLEAF  = is_leaf,
      _RULESDATA = None,
    );
    self._parser.prepare_context(context)

    # match
    processed_rules = []
    if self._struct_tags == "anyQual":
      for qname in [None] + list(quals.keys()) + ["parent"]:
        leaftag = "/".join(filter(None, [feat.type, qname]))
        leafvalue = None
        if qname:
          if qname.lower() == "parent": # because we have preprocessed parent IDS
            leafvalue = context.get("_PARENTID")
          elif qname in quals:
            leafvalue = quals[qname]
        context.tag(leaftag)
        context.update(
          force_clean = True,
          _QNAME = qname,
          _LEAFTAG = leaftag,
          _LEAFVALUE = leafvalue,
        )
        processed_rules += self._parser.process(context, ignore_unseen = True)
    elif self._struct_tags == "fullPath":
        if is_leaf:
          # processing only leaf nodes
          context.tag(fulltag)
          processed_rules += self._parser.process(context)
    else:
      raise Exception("unknown struct tags type: %s" % self._struct_tags)

    context.update_processed_rules(processed_rules)
    if processed_rules:
      self._parser.prepare_postponed(context)

    # process sub features
    if self.supporting(feat, "sub_features"):
      parent_ctx_snap = context.snap()
      ctx_parent_part = {
        "_PARENT" : feat,
        "_PARENTID" : feat.id,
        "_PARENTCTX" : parent_ctx_snap,
        "_FULLTAG" : fulltag,
      }
      for subfeat in feat.sub_features:
        context.update(ctx_parent_part, force_clean = True)
        self.process_feature(subfeat, context)

    # run postponed / actual processing on unwinding at level 1
    if depth == 1:
      self._parser.run_postponed(context)

    return

