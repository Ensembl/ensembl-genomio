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


class compare_fasta(eHive.BaseRunnable):

    def param_defaults(self):
        return {
        }

    def run(self):
        fasta1 = self.param_required("fasta1")
        fasta2 = self.param_required("fasta2")
        map_dna_path = self.param_required("seq_regions")
        output_dir = self.param_required("output_dir")
        species = self.param_required("species")
        name = self.param_required("comparison_name")
        
        map_dna = self.get_map(map_dna_path)
        seq1 = self.get_fasta(fasta1, map_dna)
        seq2 = self.get_fasta(fasta2, map_dna)
        
        comparison = self.compare_ids(seq1, seq2)
        
        output_file = output_dir + "/" + species + "_" + name + ".log"
        print("Write results in %s" % output_file)
        with open(output_file, "w") as out_fh:
            for line in comparison:
                out_fh.write(line + "\n")
        
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

    def compare_ids(self, seq1, seq2):
        comp = []
        
        # Compare number of sequences
        if len(seq1) != len(seq2):
            comp.append("WARNING: Different number of sequences: %d vs %d" % (len(seq1), len(seq2)))
        else:
            comp.append("Same number of sequences: %d" % len(seq1))
        
        # Compare all ids
       # ids1 = frozenset(seq1.keys())
       # ids2 = frozenset(seq2.keys())
       # 
       # common = frozenset.intersection(ids1, ids2)
       # comp.append("\nCommon ids: %d" % len(common))
       # 
       # diff1 = frozenset.difference(ids1, ids2)
       # if diff1:
       #     comp.append("WARNING: Ids only in 1: %d" % len(diff1))
       # diff2 = frozenset.difference(ids2, ids1)
       # 
       # if diff1:
       #     comp.append("WARNING: Ids only in 2: %d" % len(diff2))
        
        # Compare sequences
        seqs1 = {value: key for key, value in seq1.items()}
        seqs2 = {value: key for key, value in seq2.items()}
        
        common = {value for key, value in seqs1.items() if key in seqs2}
        comp.append("\nCommon sequences: %d" % len(common))

        only1 = {key: value for key, value in seqs1.items() if not key in seqs2}
        only2 = {key: value for key, value in seqs2.items() if not key in seqs1}
        
        if only1:
            total = sum([len(seq) for seq in only1.keys()])
            comp.append("WARNING: Sequences only in 1: %d (%d)" % (len(only1), total))
            only_seq1 = {name: len(seq) for seq, name in only1.items()}
            for name in sorted(only_seq1):
                comp.append("\tOnly in 1: %s (%d)" % (name, only_seq1[name]))
                
        if only2:
            total = sum([len(seq) for seq in only2.keys()])
            comp.append("WARNING: Sequences only in 2: %d (%d)" % (len(only2), total))
            only_seq2 = {name: len(seq) for seq, name in only2.items()}
            for name in sorted(only_seq2):
                comp.append("\tOnly in 2: %s (%d)" % (name, only_seq2[name]))
        
        return comp
