#!env python3

import eHive
import gzip
import io
import json
import re
import sys

from BCBio import GFF
from Bio import SeqIO
from os import path
from math import floor


class Integrity(eHive.BaseRunnable):

    def param_defaults(self):
        return {
            "ensembl_mode" : False,
        }

    def run(self):
        manifest_path = self.param_required("manifest")

        errors = []

        with open(manifest_path) as manifest_file:
            manifest = json.load(manifest_file)
            
            # Use dir name from the manifest
            for name in manifest:
                if "file" in manifest[name]:
                    file_name = manifest[name]["file"]
                    manifest[name] = path.join(path.dirname(manifest_path), file_name)
                else:
                    for f in manifest[name]:
                        if "file" in manifest[name][f]:
                            file_name = manifest[name][f]["file"]
                            manifest[name][f] = path.join(path.dirname(manifest_path), file_name)
            
            # Get content
            dna = {}
            pep = {}
            seq_regions = {}
            seq_lengths = {}
            gff = {}
            func_ann = {}
            agp_seqr = {}
            genome = {}

            if "gff3" in manifest:
                print("Got a gff")
                gff = self.get_gff3(manifest["gff3"])
            if "fasta_dna" in manifest:
                print("Got a fasta dna")
                dna, dna_errors = self.get_fasta_lengths(manifest["fasta_dna"])
                errors += dna_errors
            if "fasta_pep" in manifest:
                print("Got a fasta pep")
                pep, pep_errors = self.get_fasta_lengths(manifest["fasta_pep"])
                errors += pep_errors
            if "seq_region" in manifest:
                print("Got a seq_regions")
                seq_regions = self.get_json(manifest["seq_region"])
                seqr_lengths = {}
                seqr_seqlevel = {}
                for seq in seq_regions:
                    seq_lengths[seq["name"]] = int(seq["length"])
                    if seq["coord_system_level"] == "contig":
                        seqr_seqlevel[seq["name"]] = int(seq["length"])
            if "functional_annotation" in manifest:
                print("Got a func_anns")
                func_ann = self.get_functional_annotation(manifest['functional_annotation'])
            if "agp" in manifest:
                print("Got agp files")
                agp_seqr = self.get_agp_seq_regions(manifest['agp'])
            if "genome" in manifest:
                print("Got a genome")
                genome = self.get_json(manifest["genome"])

            # Check genome
            if genome:
                if "assembly" in genome:
                    genome_ass = genome["assembly"]
                    if "accession" in genome_ass:
                        genome_acc = genome_ass["accession"]
                        if not re.match("GCA_\d{9}(\.\d+)?", genome_acc):
                            errors += ["Genome assembly accession is wrong: '%s'" %genome_acc]

            # Check gff3
            if gff:
                if pep:
                    # We don't compare the peptide lengths because of seqedits
                    errors += self.check_lengths(pep, gff["translations"], "Fasta translations vs gff", special_diff = True)
                if func_ann:
                    errors += self.check_lengths(func_ann["genes"], gff["genes"], "Gene ids metadata vs gff", ok = "1in2")
                    errors += self.check_lengths(func_ann["translations"], gff["translations"], "Translation ids metadata vs gff", ok = "1in2")
                if seq_regions:
                    errors += self.check_seq_region_lengths(seq_lengths, gff["seq_region"], "Seq_regions metadata vs gff")

            # Check fasta dna and seq_region integrity
            if dna and seq_regions:
                errors += self.check_seq_region_lengths(seq_lengths, dna, "seq_regions json vs dna")

            # Check agp and seq_region integrity
            if agp_seqr and seq_lengths:
                errors += self.check_seq_region_lengths(seq_lengths, agp_seqr, "seq_regions json vs agps")

        if errors:
            errors_str = "\n".join(errors)
            raise Exception("Integrity test failed:\n%s" % errors_str)

    def get_fasta_lengths(self, fasta_path):
        data = {}
        non_unique = {}
        non_unique_count = 0
        empty_id_count = 0
        contains_stop_codon = 0
        for rec in SeqIO.parse(fasta_path, "fasta"):
            # Flag empty ids
            if rec.id == "":
                empty_id_count += 1
            else:
                # Flag redundant ids
                if rec.id in data:
                    non_unique[rec.id] = 1
                    non_unique_count += 1
                # Store sequence id and length
                data[rec.id] = len(rec.seq)
                if "*" in rec.seq:
                    contains_stop_codon += 1
        
        errors = []
        if empty_id_count > 0:
            errors.append("%d sequences with empty ids in %s" % (empty_id_count, fasta_path))
        if non_unique_count > 0:
            errors.append("%d non unique sequence ids in %s" % (non_unique_count, fasta_path))
        if contains_stop_codon > 0:
            errors.append("%d sequences with stop codons in %s" % (contains_stop_codon, fasta_path))
        return data, errors

    def get_json(self, json_path):
        with open(json_path) as json_file:
            return json.load(json_file)

    def get_functional_annotation(self, json_path):
        with open(json_path) as json_file:
            data = json.load(json_file)

            # Get gene ids and translation ids
            genes = {}
            translations = {}

            for item in data:
                if item["object_type"] == "gene":
                    genes[item["id"]] = 1
                elif item["object_type"] == "translation":
                    translations[item["id"]] = 1

            return { "genes" : genes, "translations" : translations }

    def get_gff3(self, gff3_path):
        if gff3_path.endswith(".gz"):
            with io.TextIOWrapper(gzip.open(gff3_path, "r")) as gff3_handle:
                return self.parse_gff3(gff3_handle)
        else:
            with open(gff3_path, "r") as gff3_handle:
                return self.parse_gff3(gff3_handle)

    def parse_gff3(self, gff3_handle):
        ensembl_mode = self.param("ensembl_mode")
        seqs = {};
        genes = {};
        peps = {};

        gff = GFF.parse(gff3_handle)
        for seq in gff:
            seqs[seq.id] = len(seq.seq)
            
            for feat in seq.features:
                if feat.type in ["gene", "ncRNA_gene", "pseudogene"]:
                    gene_id = feat.id
                    if ensembl_mode:
                        gene_id = gene_id.replace("gene:", "")
                    # Store gene length
                    genes[gene_id] = abs(feat.location.end - feat.location.start)
                    # Get CDS
                    for feat2 in feat.sub_features:
                        if feat2.type in ("mRNA", "pseudogenic_transcript"):
                            length = {}
                            for feat3 in feat2.sub_features:
                                if feat3.type == "CDS":
                                    pep_id = feat3.id
                                    if ensembl_mode:
                                        pep_id = pep_id.replace("CDS:", "")
                                    if pep_id not in length:
                                        length[pep_id] = 0
                                    length[pep_id] += abs(feat3.location.end - feat3.location.start)
                            for pep_id in length:
                                peps[pep_id] = floor(length[pep_id] / 3) - 1

        return { "seq_region": seqs, "genes": genes, "translations": peps }

    def get_agp_seq_regions(self, agp_dict):

        seqr = {}
        for agp in agp_dict:
            agp_path = agp_dict[agp]
            
            with open(agp_path, "r") as agph:
                for line in agph:
                    (asm_id, asm_start, asm_end, asm_part, typ, cmp_id, cmp_start, cmp_end, cmp_strand) = line.split("\t")
                    if typ != "W":
                        continue

                    # Assembled seq length
                    if not asm_id in seqr or seqr[asm_id] < int(asm_end):
                        seqr[asm_id] = int(asm_end)
    
                    # Composite seq length
                    if not cmp_id in seqr or seqr[cmp_id] < int(cmp_end):
                        seqr[cmp_id] = int(cmp_end)

        return seqr


    def check_ids(self, list1, list2, name):
        only1 = [];
        only2 = [];
        common = [];

        for item_id in list1:
            if item_id in list2:
                common.append(item_id)
            else:
                only1.append(item_id)
        for item_id in list2:
            if item_id not in common:
                only2.append(item_id)

        errors = []
        if common:
            print("%d common elements in %s" % (len(common), name))
        if only1:
            errors.append("%d only in first list in %s (first: %s)" % (len(only1), name, only1[0]))
        if only2:
            errors.append("%d only in second list in %s (first: %s)" % (len(only2), name, only2[0]))

        return errors

            
    def check_lengths(self, list1, list2, name, allowed_len_diff = None, ok = None, special_diff = False):
        # check list diffferences, checks if abs(values diff) < allowed_len_diff
        #  use 'ok' set to "1in2" or "2in1" if it's ok to have asymmetry 
        set1 = frozenset(list1)
        set2 = frozenset(list2)

        list1_2 = list(set1 - set2)
        list2_1 = list(set2 - set1)

        errors = []
        if len(list1_2) > 0 and (ok != "2in1"):
            errors.append("%d from the first list only for %s (i.e. %s)" % (len(list1_2), name, list1_2[0]))
        if len(list2_1) > 0 and (ok != "1in2"):
            errors.append("%d from the second list only for %s (i.e. %s)" % (len(list2_1), name, list2_1[0]))

        common_len = 0
        if allowed_len_diff is None:
            common_len = len(set1 & set2)
        else:
            # check for the sequence length difference
            diff_len_list = []
            diff_len_special_list = []
            for e in set1 & set2:
                dl12 = list1[item_id] - list2[item_id]
                if abs(dl12) <= allowed_len_diff:
                    common_len += 1
                else:
                    _dlist = diff_len_list
                    # Special case: 1 AA /BP shorter,
                    #   so assuming the stop codon is not included in the CDS (when it should be)
                    if dl12 == 1 and special_diff:
                        _dlist = diff_len_special_list
                    _dlist.append("%s: %d vs %d" % (item_id, list1[item_id], list2[item_id]))
            if diff_len_special_list:
                errors.append("%d common elements with one BP/AA length diff for %s (e.g. %s)" % (
                    len(diff_len_special_list), name, diff_len_special_list[0]
                ))
            if diff_len_list:
                errors.append("%d common elements with length diff for %s (e.g. %s)" % (
                    len(diff_len_list), name, diff_len_list[0]
                ))
        if common_len > 0:
            print("%d common elements between lists for %s" % (common_len, name), file = sys.stderr)

        return errors

    def check_seq_region_lengths(self, seqrs, feats, name):
        only_seqr = [];
        only_feat = [];

        common = [];
        diff = [];
        diff_list = [];

        for seq_id in seqrs:
            if seq_id in feats:
                # Check that feature is within the seq_region length
                if feats[seq_id] > seqrs[seq_id]:
                    diff.append(seq_id)
                    diff_list.append("%s: %d vs %d" % (seq_id, seqrs[seq_id], feats[seq_id]))
                else:
                    common.append(seq_id)
            else:
                only_seqr.append(seq_id)
        for seq_id in feats:
            if seq_id not in common and seq_id not in diff:
                only_feat.append(seq_id)

        errors = []
        if common:
            print("%d common elements in %s" % (len(common), name))
        if diff:
            errors.append("%d common elements with higher length in %s (e.g. %s)" % (len(diff), name, diff_list[0]))
        if only_seqr:
            # Not an error!
            print("%d only in seq_region list in %s (first: %s)" % (len(only_seqr), name, only_seqr[0]))
        if only_feat:
            errors.append("%d only in second list in %s (first: %s)" % (len(only_feat), name, only_feat[0]))

        return errors
