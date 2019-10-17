#!env python3

from BCBio import GFF

import eHive
import json
from Bio import SeqIO
from os import path
import gzip
import io
from math import floor

class Integrity(eHive.BaseRunnable):

    def param_default(self):
        return {
        }

    def run(self):
        manifest_path = self.param_required("manifest")
        #print(manifest_path)

        errors = []

        with open(manifest_path) as manifest_file:
            manifest = json.load(manifest_file)
            
            # Use dir name from the manifest
            for f in manifest:
                manifest[f] = path.join(path.dirname(manifest_path), manifest[f])
            
            # Get content
            dna = {}
            pep = {}
            seq_regions = {}
            seq_lengths = {}
            gff = {}

            if "gff3" in manifest:
                gff = self.get_gff3(manifest["gff3"])
            if "fasta_dna" in manifest:
                dna = self.get_fasta_lengths(manifest["fasta_dna"])
            if "fasta_pep" in manifest:
                pep = self.get_fasta_lengths(manifest["fasta_pep"])
            if "metadata_seq_region" in manifest:
                seq_regions = self.get_json(manifest["metadata_seq_region"])
                seqr_lengths = {}
                for seq in seq_regions:
                    seq_lengths[seq["name"]] = seq["length"]

            # Check gff3
            if gff:
                if pep:
                    # We don't compare the peptide lengths because of seqedits
                    errors += self.check_lengths(pep, gff["pep"], "fasta peptide vs gff CDS")
                if seq_regions:
                    errors += self.check_lengths(seq_lengths, gff["seq_region"], "seq_regions json vs gff seq_regions", 0)

            # Check fasta dna and seq_region integrity
            if dna and seq_regions:
                # We don't tolerate any length difference
                errors += self.check_lengths(dna, seq_lengths, "fasta dna vs seq_regions json", 0)

        if errors:
            raise Exception("Integrity test failed: " + str(errors))

    def get_fasta_lengths(self, fasta_path):
        data = {}
        for rec in SeqIO.parse(fasta_path, "fasta"):
            data[rec.id] = len(rec)
        return data

    def get_json(self, json_path):
        with open(json_path) as json_file:
            return json.load(json_file)

    def get_gff3(self, gff3_path):

        seqs = {};
        genes = {};
        peps = {};

        with io.TextIOWrapper(gzip.open(gff3_path, "r")) as gff3_handle:
            gff = GFF.parse(gff3_handle)
            for seq in gff:
                seqs[seq.id] = len(seq.seq)
                
                for feat in seq.features:
                    if feat.type == "gene":
                        genes[feat.id] = abs(feat.location.end - feat.location.start)
                        # Get CDS
                        for feat2 in feat.sub_features:
                            if feat2.type == "mRNA":
                                length = {}
                                for feat3 in feat2.sub_features:
                                    if feat3.type == "CDS":
                                        pep_id = feat3.id.replace("CDS:", "")
                                        if pep_id not in length:
                                            length[pep_id] = 0
                                        length[pep_id] += abs(feat3.location.end - feat3.location.start)
                                for pep_id in length:
                                    peps[pep_id] = floor(length[pep_id] / 3) - 1

        return { "seq_region": seqs, "gene": genes, "pep": peps }

            
    def check_lengths(self, list1, list2, name, allowed_diff = None):
        only1 = [];
        only2 = [];

        common = [];
        diff = [];
        diff_list = [];
        diff1 = [];
        diff1_list = [];

        for item_id in list1:
            if item_id in list2:
                # Check that the difference of length is within a threshold
                if allowed_diff is not None and abs(list1[item_id] - list2[item_id]) > allowed_diff:
                    # Special case: shorter by one, so assuming the stop codon is not included in the CDS (when it should be)
                    if list1[item_id] - list2[item_id] == 1:
                        diff1.append(item_id)
                        diff1_list.append("%s: %d vs %d" % (item_id, list1[item_id], list2[item_id]))
                    else:
                        diff.append(item_id)
                        diff_list.append("%s: %d vs %d" % (item_id, list1[item_id], list2[item_id]))
                else:
                    common.append(item_id)
            else:
                only1.append(item_id)
        for item_id in list2:
            if item_id not in common and item_id not in diff and item_id not in diff1:
                only2.append(item_id)

        errors = []
        if common:
            print("%d common elements in %s" % (len(common), name))
        if diff1:
            errors.append("%d common elements with one shorter in %s (e.g. %s)" % (len(diff1), name, diff1_list[0]))
        if diff:
            errors.append("%d common elements with different length in %s (e.g. %s)" % (len(diff), name, diff_list[0]))
        if only1:
            errors.append("%d only in first list in %s (first: %s)" % (len(only1), name, only1[0]))
        if only2:
            errors.append("%d only in second list in %s (first: %s)" % (len(only2), name, only2[0]))

        return errors

    def check_gff_seq_regions(gff, seqr):
        return

