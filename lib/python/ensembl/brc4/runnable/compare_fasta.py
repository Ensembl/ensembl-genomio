#!env python3

import eHive
import gzip
import json
import io
import re
import sys
from functools import partial

from Bio import SeqIO
from os import path
from ensembl.brc4.runnable.parser import Parser


class compare_fasta(eHive.BaseRunnable):

    def param_defaults(self):
        return {
        }

    def run(self):
        report = self.param_required("report")
        fasta1 = self.param_required("fasta1")
        fasta2 = self.param_required("fasta2")
        map_dna_path = self.param_required("seq_regions")
        output_dir = self.param_required("output_dir")
        species = self.param_required("species")
        name = self.param_required("comparison_name")

        map_dna = self.get_map(map_dna_path)
        seq1 = self.get_fasta(fasta1, map_dna)
        seq2 = self.get_fasta(fasta2, map_dna)

        (stats, diffs, seq_map) = self.compare_seqs(seq1, seq2)

        # Print mapping to a file (add report data)
        map_file = output_dir + "/" + species + "_" + name + ".map"
        self.print_map(seq_map, map_file, report)

        # Print full list of results in a file
        output_file = output_dir + "/" + species + "_" + name + ".log"
        print("Write results in %s" % output_file)
        with open(output_file, "w") as out_fh:
            for line in diffs:
                out_fh.write(line + "\n")

        # Print the stats separately
        out = {
            "species": species,
            "stats": stats
        }
        self.dataflow(out, 2)

    def print_map(self, seq_map, map_file, report_file):

        report_parser = Parser()
        report_seq = report_parser.get_report_regions(report_file)
        report = self.add_report_to_map(seq_map, report_seq)

        print("Write map in %s" % map_file)
        with open(map_file, "w") as out_fh:
            out_fh.write(json.dumps(report, sort_keys=True, indent=4))

    def add_report_to_map(self, seq_map, report_seq):

        accession_version = r'\.\d+$'
        report = []
        for old_name, insdc_name in seq_map.items():
            if insdc_name not in report_seq:
                raise Exception("No INSDC %s found in report" % insdc_name)
            else:
                seqr = report_seq[insdc_name]
                seqr["name"] = old_name
                seqr["EBI_seq_region_name"] = old_name
                brc4_name = insdc_name
                brc4_name = re.sub(accession_version, '', brc4_name)
                seqr["BRC4_seq_region_name"] = brc4_name
                syns = [{
                    "source": "INSDC",
                    "name": insdc_name
                }]
                seqr["synonyms"] = syns
                report.append(seqr)
        return report

    def get_map(self, map_path):

        print("Read file %s" % map_path)
        data = self.get_json(map_path)

        map_dna = {}

        for seqr in data:
            if "synonyms" in seqr:
                for syn in seqr["synonyms"]:
                    if syn["name"] == "INSDC":
                        map_dna[name] = syn["value"]

        return map_dna

    def get_json(self, json_path):
        with open(json_path) as json_file:
            return json.load(json_file)

    def get_fasta(self, fasta_path, map_dna):

        print("Read file %s" % fasta_path)
        sequences = {}
        _open = partial(gzip.open, mode='rt') if fasta_path.endswith(
            '.gz') else open
        with _open(fasta_path) as fasta_fh:
            for rec in SeqIO.parse(fasta_fh, "fasta"):
                name = rec.id
                if name in map_dna:
                    name = map_dna[name]
                sequences[name] = re.sub(r"[^CGTA]", "N", str(rec.seq.upper()))
        return sequences

    def compare_seqs(self, seq1, seq2):
        comp = []
        diff = abs(len(seq1) - len(seq2))
        stats = {
            "seq_count_1": len(seq1),
            "seq_count_2": len(seq2),
            "num_diff_seq": diff,
            "common": 0,
            "only1": 0,
            "only2": 0,
            "max_only1": 0,
            "max_only2": 0,
            "only1_200": 0,
            "only1_1000": 0,
            "only2_200": 0,
            "only2_1000": 0,
            "other_locations": 0,
            "summary": None,
            "organellar_summary": None
        }

        value = None  # variable used for summary
        org_value = None  # variable used for organellar_summary

        # Compare number of sequences
        if len(seq1) != len(seq2):
            comp.append("WARNING: Different number of sequences: %d vs %d" % (
                len(seq1), len(seq2)))
            value = "mismatch"
            org_value = "unkown"
        else:
            comp.append("Same number of sequences: %d" % len(seq1))
            value = "identical"

        # Compare sequences
        seqs1 = {seq: name for name, seq in seq1.items()}
        seqs2 = {seq: name for name, seq in seq2.items()}

        common = {seqs2[seq]: name for seq,
                  name in seqs1.items() if seq in seqs2}
        only1 = {seq: name for seq, name in seqs1.items() if not seq in seqs2}
        only2 = {seq: name for seq, name in seqs2.items() if not seq in seqs1}

        # comparing the organellar sequences
        report = self.param_required("report")
        report_parser = Parser()
        report_seq = report_parser.get_report_regions(report)
        report = self.add_report_to_map(common, report_seq)
        map_dna_path = self.param_required("seq_regions")
        data = self.get_json(map_dna_path)
        org = []
        org1 = []
        org2 = []
        org3 = []

        # Gathering data from the INSDC report file and storing it into a list
        for k, v in report_seq.items():
            if "location" in v:
                if v['location'] != "chromosome" and v['location'] != "nuclear_chromosome" and v['location'] != "linkage_group":
                    x = v['location']
                    org1.append(x)
                    org.append(k)
            else:
                pass

        # Gathering data from Seq_json file and storing it into a list
        for j in data:
            for k, v in j.items():
                if "location" in k:
                    if v != "chromosome" and v != "nuclear_chromosome" and v != "linkage_group":
                        y = j['EBI_seq_region_name']
                        org3.append(v)
                        org2.append(y)
                else:
                    pass

        # checking for INSDC names for organellar sequences found in core
        for i in range(0, len(org2)):
            for p, l in report_seq.items():
                if l['name'] == org2[i]:
                    org2[i] = p
        for index, i in enumerate(org2):
            if i not in org:
                org.append(i)
                org1.append(org3[index])
            else:
                print("you have duplicates")

        # checking for multiple entries of organellar seq
        myList = [i.split('.')[0] for i in org]
        myList1 = [j[:-1] for j in myList]  # similar accession
        a = list(set(myList1))

        # comparing organellar sequences with common, only1 and only2
        count = 0
        for index, i in enumerate(org):
            if i == 'na':
                comp.append("MISSING acession in the report (na)")
            else:
                if i in common.values():
                    count = count+1
                    comp.append("%s (both) in location: %s" % (i, org1[index]))
                    if value == "identical" and count > 0:
                        org_value = "organell_added"
                elif i in only1.values():
                    count = count+1
                    comp.append("%s (only1) in  location: %s" %
                                (i, org1[index]))
                    org_value = "unkown_with_organelle"
                else:
                    count = count+1
                    comp.append("%s (only2) in location: %s" %
                                (i, org1[index]))
                    org_value = "unkown_with_organelle"

        # if the mistmatch is due to added organell
        if len(seq1) > len(seq2):
            greater_len = len(seq1)
        else:
            greater_len = len(seq2)
        diff_common = greater_len - len(common)

        if diff != 0:
            if diff == count and diff_common == count:
                org_value = "organell_added"

        # checking if multiple entries of organellar chromosomes are present
        if len(myList1) != len(a):
            org_value = "WARNING:Multiple_entry"

        # updating the stats
        stats["common"] = len(common)
        stats["only1"] = len(only1)
        stats["only2"] = len(only2)
        stats["other_locations"] = count
        stats["summary"] = value
        stats["organellar_summary"] = org_value

        if len(only1) > 0 or len(only2) > 0:
            comp.append("\nCommon sequences: %d" % len(common))

            if only1:
                stats["max_only1"] = max(len(seq) for seq in only1.keys())

                # Only list sequences where the length is > 200
                mini = {seq: name for seq, name in only1.items()
                        if len(seq) <= 200}
                maxi = {seq: name for seq, name in only1.items()
                        if len(seq) > 200}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 1 with length <= 200" % len(mini))
                    only1 = maxi

            if only1:
                # Only list sequences where the length is > 1000
                mini = {seq: name for seq, name in only1.items()
                        if len(seq) <= 1000}
                maxi = {seq: name for seq, name in only1.items()
                        if len(seq) > 1000}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 1 with length <= 1000" % len(mini))
                    only1 = maxi

            if only1:
                total = sum([len(seq) for seq in only1.keys()])
                comp.append("WARNING: Sequences only in 1: %d (%d)" %
                            (len(only1), total))
                only_seq1 = {name: len(seq) for seq, name in only1.items()}
                for name, length in sorted(only_seq1.items(), key=lambda x: x[1]):
                    comp.append("\tOnly in 1: %s (%d)" % (name, length))

            if only2:
                stats["max_only2"] = max(len(seq) for seq in only2.keys())

                # Only list sequences where the length is > 200
                mini = {seq: name for seq, name in only2.items()
                        if len(seq) <= 200}
                maxi = {seq: name for seq, name in only2.items()
                        if len(seq) > 200}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 2 with length <= 200" % len(mini))
                    only2 = maxi

            if only2:
                # Only list sequences where the length is > 1000
                mini = {seq: name for seq, name in only2.items()
                        if len(seq) <= 1000}
                maxi = {seq: name for seq, name in only2.items()
                        if len(seq) > 1000}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 2 with length <= 1000" % len(mini))
                    only2 = maxi

            if only2:
                total = sum([len(seq) for seq in only2.keys()])
                comp.append("WARNING: Sequences only in 2: %d (%d)" %
                            (len(only2), total))
                only_seq2 = {name: len(seq) for seq, name in only2.items()}
                for name, length in sorted(only_seq2.items(), key=lambda x: x[1]):
                    comp.append("\tOnly in 2: %s (%d)" % (name, length))

        return (stats, comp, common)
