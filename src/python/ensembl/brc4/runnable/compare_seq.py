import argparse
import logging
import hashlib
from pathlib import Path
import re
from os import PathLike
from typing import List, Any, Tuple

from Bio import SeqIO

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


class SeqGroup:
    def __init__(self, sequence, identifier=None) -> None:
        self.sequence = sequence
        self.length = len(self.sequence)
        self.ids = []
        if identifier:
            self.add_id(identifier)
        self.count = len(self.ids)

    def __str__(self) -> str:
        return ", ".join(self.ids)

    def add_id(self, identifier) -> None:
        self.ids.append(identifier)
        self.count = len(self.ids)

class CompareFasta:
    def __init__(self, fasta1_path: str, fasta2_path: str, output_dir) -> None:
        self.fasta1 = fasta1_path
        self.fasta2 = fasta2_path
        self.output_dir = output_dir
        self.common = {}
        self.comp = []
        self.compare_seqs(fasta1_path, fasta2_path)

    def read_fasta(self, fasta_path) -> dict:
        print(f"Read file {fasta_path}")
        sequences = {}
        with open_gz_file(fasta_path) as fasta_fh:
            for rec in SeqIO.parse(fasta_fh, "fasta"):
                name = rec.id
                #if name in map_dna:      ---> confirm if ensembl ever uses different name for loading?
                #    name = map_dna[name]
                sequences[name] = re.sub(r"[^CGTA]", "N", str(rec.seq.upper()))
        return sequences
    
    def build_seq_dict(self, seqs: dict) -> dict:
        """Build a seq dict taking duplicates into account"""

        seqs_dict = dict()
        for name, seq in seqs.items():
            if seq in seqs_dict:
                seqs_dict[seq].add_id(name)
            else:
                seqs_dict[seq] = SeqGroup(seq, name)
        return seqs_dict
    
    def find_common_groups(self, seqs1: dict, seqs2: dict) -> Tuple[dict, List[Any]]:
        common = {}
        for seq1, group1 in seqs1.items():
            if seq1 in seqs2:
                group2 = seqs2[seq1]
                # Check that the 2 groups have the same number of sequences
                if group1.count == group2.count:
                    if group1.count == 1:
                        common[group1.ids[0]] = group2.ids[0]
                    else:
                        self.comp.append(f"Matched 2 identical groups of sequences: {group1} and {group2}")
                        possible_id2 = " OR ".join(group2.ids)
                        for id1 in group1.ids:
                            common[id1] = possible_id2
                else:
                   self.comp.append(
                        f"Matched 2 different groups of sequences ({group1.count} vs {group2.count}): {group1} and {group2}"
                    )

        return common,self.comp
    
    def compare_seqs(self, fasta1, fasta2):
        """Todo"""
        seq1 = self.read_fasta(fasta1)
        seq2 = self.read_fasta(fasta2)

        # Compare sequences
        seq1_dict = self.build_seq_dict(seq1)
        seq2_dict = self.build_seq_dict(seq2)

        # Compare number of sequences
        if len(seq1) != len(seq2):
           self.comp.append(f"WARNING: Different number of sequences: {len(seq1)} vs {len(seq2)}")

        
        common, group_comp = self.find_common_groups(seq1_dict, seq2_dict)
        self.comp += group_comp
        
        # Sequences that are not common
        only1 = {seq: group for seq, group in seq1_dict.items() if not seq in seq2_dict}

        only2 = {seq: group for seq, group in seq2_dict.items() if not seq in seq1_dict}
        
        if only1:
            for ids in only1.values():
                self.comp.append(f"Only in fasta1: {ids}\n")
        if only2:
            for ids in only2.values():
                self.comp.append(f"Only in fasta2: {ids}\n")
        if common:
            for common_ids in common.items():
                self.comp.append(f"Common ids: {common_ids}")
        
        self.check_for_N(only1, only2)
        # Print full list of results in a file
        output_file = Path.joinpath(self.output_dir, "compare.log")
        print(f"Write results in {output_file}")
        with open(output_file, "w") as out_fh:
            for line in self.comp:
                out_fh.write(line + "\n")
    
    def check_for_N(self, only1, only2):
        names_length = {}
        # sequences which have extra N at the end
        if only1 and only2:
            for seq_1, name1 in only1.items():
                len1 = len(seq_1)
                seq1_N = seq_1.count("N")

                for seq_2, name2 in only2.items():
                    len2 = len(seq_2)
                    seq2_N = seq_2.count("N")

                    if seq_2[:len1] == seq_1:
                        ignored_seq = seq_2[len1:]
                        N = ignored_seq.count("N")

                        if len(ignored_seq) == N:
                           self.comp.append(f"Please check extra Ns added in core in {name1} and {name2}")
                        else:
                           self.comp.append(
                                f"ALERT INSERTIONS at the end or diff assembly level {name1} and {name2}"
                            )

                    elif len1 == len2:
                        if seq2_N > seq1_N:
                           self.comp.append(f"Core has more Ns, check {name1} and {name2}")
                        elif seq1_N > seq2_N:
                           self.comp.append(f"INSDC has more Ns, check {name1} and {name2}")
                        else:
                            names_length[name1] = name2

        if names_length:
            length = len(names_length)
            self.comp.append(f"{length} sequences have the same length")
            for insdc, core in names_length.items():
               self.comp.append(f"INSDC: {insdc} and coredb : {core}")

def parse_args(arg_list: list[str] | None) -> argparse.Namespace:
    """TODO
    Args:
    arg_list: TODO
    """
    parser = ArgumentParser(description="Compare sequences between two genomes")
    # Add filter arguments
    parser.add_argument("--fasta1_path", required=True, help="Path to INSDC fasta file")
    parser.add_argument("--fasta2_path", required=True, help="Path to user fasta file")
    parser.add_argument("--output_dir", default=Path.cwd(), help="Folder to store the report files")
                        
    # Add flags
    parser.add_argument(
        "--compare_seq_region",
            action="store_true",
        help="Enable compare seq_region mode, i.e. use seq_region for seqence comparison",
    )
    return parser.parse_args(arg_list)

def main(arg_list: list[str] | None = None) -> None:
    """Main script entry-point.
    Args:
    arg_list: TODO
     """
    args = parse_args(arg_list)

    CompareFasta(args.fasta1_path, args.fasta2_path, args.output_dir)

if __name__ == "__main__":
    main()