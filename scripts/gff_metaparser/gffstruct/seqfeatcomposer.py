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


import sys

from BCBio import GFF

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


class SeqFeatComposer:
    # storing/building gff structures
    def __init__(self):
        self._topf = set()
        self._processed = {}
        self._supported_fields = dict()

    def gff_write(
        self, out_rec_features, contig_id=None, contig_len=None, out_file=None, global_context=None
    ):
        if out_file is None:
            return
        if out_rec_features:
            self.update_stashed(out_rec_features, global_context)
            out_rec = SeqRecord(Seq(None, length=contig_len), id=contig_id)
            out_rec.features = out_rec_features
            GFF.write([out_rec], out_file)

    def update_stashed(self, features, global_context):
        if not features or not global_context:
            return
        for feat in features:
            quals = self.supporting(feat, "qualifiers") and feat.qualifiers or {}
            del_list = []
            for q in quals.keys():
                val = quals[q]
                if type(val) == dict and "_FROM_STASH" in val:
                    val = val["_FROM_STASH"]
                    val = global_context.get(val.get("obj_tag"), val.get("id"), val.get("path"))
                    if val is None:
                        del_list.append(q)
                    else:
                        quals[q] = val
            for q in del_list:
                del quals[q]
            if self.supporting(feat, "sub_features"):
                self.update_stashed(feat.sub_features, global_context)
        return

    def supporting(self, obj, field):
        tp = type(obj)
        if tp not in self._supported_fields:
            self._supported_fields[tp] = frozenset(obj.__dir__())
        return field in self._supported_fields[tp]

    def cid(self, ctx):
        p_type, p_id, p_ctx = "", "", ctx.get("_PARENTCTX")
        if p_ctx:
            p_type, p_id = map(lambda k: p_ctx.get(k, ""), ["_TYPE", "_ID"])
        _type, _start, _end, _strand, _id = map(
            lambda k: ctx.get(k, ""), ["_TYPE", "_START", "_END", "_STRAND", "_ID"]
        )
        return (p_type, p_id, _type, _start, _end, _strand, _id)

    def processed_add(self, ctx, used_quals=None, kid=None):
        cid = self.cid(ctx)
        feat = self._processed.get(cid)

        # flattern and exclude "parent"s
        used_quals_flat = {v[0]: v[1] for k, v in (used_quals or {}).items() if k != "parent" and k[0] != "_"}

        if not feat:
            location = ctx["_LOCATION"]  # die if no location
            _type = ctx["_TYPE"]  # should be already changed by rules
            is_leaf = ctx["_ISLEAF"]
            strand = ctx.get("_STRAND")
            phase = ctx.get("_PHASE")
            _id = ctx.get("_ID")
            source = ctx.get("_SRC")
            # quals
            quals = used_quals_flat or {}
            if not is_leaf and "ID" not in quals:
                quals["ID"] = _id
            if phase is not None and "phase" not in quals:
                quals["phase"] = phase
            # create feat object
            quals = self.sort_quals(quals)
            location.strand = strand
            obj = SeqFeature(
                location, type=_type, qualifiers=quals
            )  # NB referebce to quals / allows post assign modifications
            obj.source = source
            # fill processed
            feat = {
                "cid": cid,
                "obj": obj,
                "quals": quals,
                "kids": dict(),
                "is_leaf": is_leaf,
                "source": source,
            }
            self._processed[cid] = feat
        else:
            # fill quals
            used_quals_flat = self.sort_quals(used_quals_flat)
            # assert( id(feat["quals"]) = id(feat["obj"].qualifiers) )
            # print(cid, hex(id(feat["quals"])), hex(id(feat["obj"].qualifiers)), file=sys.stderr)
            feat["quals"].update(used_quals_flat)

        if kid:
            kcid = kid["cid"]
            if kcid not in feat["kids"]:
                feat["kids"][kcid] = len(feat["kids"])

        return feat

    def sort_quals(self, quals):
        if not quals:
            return quals
        filtered_id = filter(lambda x: x[0] == "ID", quals.items())
        filtered = filter(lambda x: x[0][0] != "_" and x[0] != "ID", quals.items())
        return dict(list(filtered_id) + list(sorted(filtered, key=lambda x: x[0])))

    def compose(self, context, out=None):
        if out is None:
            return

        # clear
        self._topf = set()
        self._processed = {}

        # go from used leaves to top
        for ctx in context.used_leaves():
            # print("compose ctx type", ctx.get("_TYPE"), file =sys.stderr)
            rules_data = ctx.get("_RULESDATA")
            if rules_data and "_ALL" in rules_data and "USEDQUALS" in rules_data["_ALL"]:
                used_quals = rules_data["_ALL"]["USEDQUALS"]

            # TODO: rewrite in a sane manner. see fix_actions
            feat = self.processed_add(ctx, used_quals=used_quals)
            parent_ctx = ctx.get("_PARENTCTX")

            while parent_ctx:
                rules_data = parent_ctx.get("_RULESDATA")
                if rules_data and "_ALL" in rules_data and "USEDQUALS" in rules_data["_ALL"]:
                    used_quals = rules_data["_ALL"]["USEDQUALS"]

                feat = self.processed_add(parent_ctx, used_quals=used_quals, kid=feat)
                parent_ctx = parent_ctx.get("_PARENTCTX")

            self._topf.add(feat["cid"])

        # update subfeatures
        for cid, feat in self._processed.items():
            if feat["kids"]:
                feat["obj"].sub_features = [
                    self._processed[cid]["obj"]
                    for cid, _ in sorted(feat["kids"].items(), key=lambda k: -k[1])
                ]
            else:
                feat["obj"].sub_features = []
        # fill out
        for cid in sorted(self._topf):
            # print(cid, file = sys.stderr)
            # print(self._processed[cid], file = sys.stderr)
            # pass
            out.append(self._processed[cid]["obj"])
