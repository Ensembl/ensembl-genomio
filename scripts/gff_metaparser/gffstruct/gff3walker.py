# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import functools
import re
import sys

from BCBio import GFF

# locals
from .seqfeatcomposer import SeqFeatComposer
from .utils import UpdatingLen
from .walkcontext import WalkContext


class GFF3Walker:
    _valid_structure_tags = frozenset("fullPath anyQual".split())  # /gene/mrna/exon vs exon/id, exon

    def __init__(self, parser, in_file, structure_tags="fullPath", global_ctx=None, norm_id=None):
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

    def norm_id(self, id, type=None, context=None):
        if self._norm_id is None:
            return id
        return self._norm_id(id, type=type, context=context)

    def norm_and_update_id(self, feat, quals, context):
        _id, _type = feat.id, feat.type
        if not _id:
            return _id
        if _type:
            _type = _type.lower()
        _raw_id = _id
        # can norm
        if self._norm_id and self._norm_id.supporting(feat.type):
            _id = self.norm_id(_id, _type, context)  # or cotext.data???
        # fix duplicates
        _id = context.fix_id_sfx(_id, _type)
        # update
        if _raw_id != _id:
            feat.id = _id
            if quals and quals.get("ID"):
                quals["ID"][0] = _id
        return _id

    def supporting(self, obj, field):
        tp = type(obj)
        if tp not in self._supported_fields:
            self._supported_fields[tp] = frozenset(obj.__dir__())
        return field in self._supported_fields[tp]

    def walk(self, out_file=None, seq_len_dict=None, contig_length_extension=True):
        gff = GFF.parse(self._in_file)
        sf_composer = SeqFeatComposer()  #  ?? init composer once, get on __init__???
        for contig in gff:
            # get len from gff or fna file, try to infere if not available
            ctg_len_pre = len(contig)
            if not contig_length_extension:
                ctg_len_pre = None
            ctg_len = UpdatingLen(ctg_len_pre, force_update=contig_length_extension)
            ctg_len.update(seq_len_dict(contig.id), stop_on_success=True)
            # feature processing
            out_rec_features = None
            if out_file:
                out_rec_features = []
            for cnt, topft in enumerate(contig.features):
                # iterate through all of the features
                context = WalkContext(global_context=self.global_context, ctg_len_inferer=ctg_len)
                context.update("_SEQID", contig.id)
                self.process_feature(topft, context)
                # compose gff
                sf_composer.compose(context, out_rec_features)
            # dump gff
            sf_composer.gff_write(
                out_rec_features,
                contig_id=contig.id,
                contig_len=len(ctg_len),
                out_file=out_file,
                global_context=context.global_context,
            )

    def process_feature(self, feat, context):
        # store stats in the context
        loc = feat.location
        quals = feat.qualifiers
        phase = quals.get("phase")
        fulltag_parts = list(filter(None, [context.get("_FULLTAG"), feat.type]))
        fulltag = "/".join(fulltag_parts)
        depth = fulltag.count("/") + 1
        source = self.supporting(feat, "source") and feat.source or None
        is_leaf = (
            not self.supporting(feat, "sub_features") or not feat.sub_features
        )  # fullPath processed and ID can be ommited only for leaves

        if depth == 1:
            # set top feature
            context.top(feat)

        context.update(
            force_clean=True,
            _SRC=source,
            _TYPE=feat.type,
            _START=loc.start,
            _END=loc.end,
            _STRAND=feat.location.strand,
            _LOCATION=feat.location,
            _PHASE=phase,
            _FULLTAG=fulltag,
            _DEPTH=depth,
            _ISLEAF=is_leaf,
            _RULESDATA=None,
        )

        # norm
        new_id = self.norm_and_update_id(feat, quals, context)

        context.update(
            force_clean=True,
            _ID=new_id,
            _QUALS=quals,
        )

        self._parser.prepare_context(context)

        # update contig length
        context.update_ctg_len(loc.end)

        # match
        processed_rules = []
        if self._struct_tags == "anyQual":
            allowed_add_ons = ["parent", "_START", "_END"]
            for qname in [None] + list(quals.keys()) + allowed_add_ons:
                leaftag = "/".join(filter(None, [feat.type, qname]))
                leafvalue = None
                if qname:
                    if qname.lower() == "parent":  # because we have preprocessed parent IDS
                        leafvalue = context.get("_PARENTID")
                    elif qname in allowed_add_ons:
                        leafvalue = context.get(qname)
                    elif qname in quals:
                        leafvalue = quals[qname]
                context.tag(leaftag)
                context.update(
                    force_clean=True,
                    _QNAME=qname,
                    _LEAFTAG=leaftag,
                    _LEAFVALUE=leafvalue,
                )
                processed_rules += self._parser.process(context, ignore_unseen=True)
                # print(fulltag, leaftag, "processed", len(processed_rules), file=sys.stderr)
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
                "_PARENT": feat,
                "_PARENTID": feat.id,
                "_PARENTCTX": parent_ctx_snap,
                "_FULLTAG": fulltag,
            }
            for subfeat in feat.sub_features:
                context.update(ctx_parent_part, force_clean=True)
                self.process_feature(subfeat, context)

        # run postponed / actual processing on unwinding at level 1
        if depth == 1:
            self._parser.run_postponed(context)

        return
