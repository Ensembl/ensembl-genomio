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
                "species" : species,
                "stats" : stats
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
        _open = partial(gzip.open, mode='rt') if fasta_path.endswith('.gz') else open
        with _open(fasta_path) as fasta_fh:
            for rec in SeqIO.parse(fasta_fh, "fasta"):
                name = rec.id
                if name in map_dna:
                    name = map_dna[name]
                sequences[name] = re.sub(r"[^CGTA]", "N", str(rec.seq.upper()))
        return sequences

    def compare_seqs(self, seq1, seq2):
        comp = []
        stats = {
                "seq_count_1" : len(seq1),
                "seq_count_2" : len(seq2),
                "diff_length" : abs(len(seq1) - len(seq2)),
                "common" : 0,
                "only1" : 0,
                "only2" : 0,
                "max_only1" : 0,
                "max_only2" : 0,
                "only1_200" : 0,
                "only1_1000" : 0,
                "only2_200" : 0,
                "only2_1000" : 0,
        }
        
        # Compare number of sequences
        if len(seq1) != len(seq2):
            comp.append("WARNING: Different number of sequences: %d vs %d" % (len(seq1), len(seq2)))
        else:
            comp.append("Same number of sequences: %d" % len(seq1))
        
        # Compare sequences
        seqs1 = {seq: name for name, seq in seq1.items()}
        seqs2 = {seq: name for name, seq in seq2.items()}
        
        common = {seqs2[seq]: name for seq, name in seqs1.items() if seq in seqs2}
        only1 = {seq: name for seq, name in seqs1.items() if not seq in seqs2}
        only2 = {seq: name for seq, name in seqs2.items() if not seq in seqs1}

        stats["common"] = len(common)
        stats["only1"] = len(only1)
        stats["only2"] = len(only2)
        
        if len(only1) > 0 or len(only2) > 0:
            comp.append("\nCommon sequences: %d" % len(common))

            if only1:
                stats["max_only1"] = max(len(seq) for seq in only1.keys())
                
                # Only list sequences where the length is > 200
                mini = {seq: name for seq, name in only1.items() if len(seq) <= 200}
                maxi = {seq: name for seq, name in only1.items() if len(seq) > 200}
                
                if mini and len(mini) > 3000:
                    comp.append("WARNING: Ignoring %d sequences from 1 with length <= 200" % len(mini))
                    only1 = maxi
                
            if only1:
                # Only list sequences where the length is > 1000
                mini = {seq: name for seq, name in only1.items() if len(seq) <= 1000}
                maxi = {seq: name for seq, name in only1.items() if len(seq) > 1000}

                if mini and len(mini) > 3000:
                    comp.append("WARNING: Ignoring %d sequences from 1 with length <= 1000" % len(mini))
                    only1 = maxi
            
            if only1:
                total = sum([len(seq) for seq in only1.keys()])
                comp.append("WARNING: Sequences only in 1: %d (%d)" % (len(only1), total))
                only_seq1 = {name: len(seq) for seq, name in only1.items()}
                for name, length in sorted(only_seq1.items(), key=lambda x: x[1]):
                    comp.append("\tOnly in 1: %s (%d)" % (name, length))
                    
            if only2:
                stats["max_only2"] = max(len(seq) for seq in only2.keys())

                # Only list sequences where the length is > 200
                mini = {seq: name for seq, name in only2.items() if len(seq) <= 200}
                maxi = {seq: name for seq, name in only2.items() if len(seq) > 200}
                
                if mini and len(mini) > 3000:
                    comp.append("WARNING: Ignoring %d sequences from 2 with length <= 200" % len(mini))
                    only2 = maxi
                    
            if only2:
                # Only list sequences where the length is > 1000
                mini = {seq: name for seq, name in only2.items() if len(seq) <= 1000}
                maxi = {seq: name for seq, name in only2.items() if len(seq) > 1000}
                
                if mini and len(mini) > 3000:
                    comp.append("WARNING: Ignoring %d sequences from 2 with length <= 1000" % len(mini))
                    only2 = maxi

            if only2:
                total = sum([len(seq) for seq in only2.keys()])
                comp.append("WARNING: Sequences only in 2: %d (%d)" % (len(only2), total))
                only_seq2 = {name: len(seq) for seq, name in only2.items()}
                for name, length in sorted(only_seq2.items(), key=lambda x: x[1]):
                    comp.append("\tOnly in 2: %s (%d)" % (name, length))
        
        return (stats, comp, common)
