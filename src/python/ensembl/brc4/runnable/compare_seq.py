import argparse
from pathlib import Path
import re
from typing import Any, List, Tuple

from Bio import SeqIO

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser


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
    def __init__(self, fasta1_path: Path, fasta2_path: Path, output_dir: str) -> None:
        self.fasta1 = Path(fasta1_path)
        self.fasta2 = Path(fasta2_path)
        self.output_dir = output_dir
        self.comp = []
        self.compare_seqs(fasta1_path, fasta2_path)

    def read_fasta(self, fasta_path: Path) -> dict:
        """
        Reads a FASTA file and returns a dictionary mapping sequence IDs to sequences.

        Args:
            fasta_path (Path): Path to the FASTA file. Supports gzipped files.

        Returns:
            dict: A dictionary where keys are sequence IDs and values are sequences
              with all non-CGTA characters replaced by 'N'.
        """
        print(f"Read file {fasta_path}")
        sequences = {}
        with open_gz_file(fasta_path) as fasta_fh:
            for rec in SeqIO.parse(fasta_fh, "fasta"):
                name = rec.id
                sequences[name] = re.sub(r"[^CGTA]", "N", str(rec.seq.upper()))
        return sequences
    
    def build_seq_dict(self, seqs: dict) -> dict:
        """
        Builds a dictionary of unique sequences and their associated IDs, accounting for duplicates.

        Args:
            seqs (dict): A dictionary where keys are sequence IDs and values are sequences.

        Returns:
            dict: A dictionary where keys are unique sequences and values are `SeqGroup` objects
                that group sequence IDs sharing the same sequence.
        """
        seqs_dict = dict()
        for name, seq in seqs.items():
            if seq in seqs_dict:
                seqs_dict[seq].add_id(name)
            else:
                seqs_dict[seq] = SeqGroup(seq, name)
        return seqs_dict
    
    def find_common_groups(self, seq1_dict: dict, seq2_dict: dict) -> Tuple[dict, List[Any]]:
        """
        Find common sequences between two dictionaries and group them.

        Args:
            seq1_dict (dict): Dictionary of sequences from the first dataset.
            seq2_dict (dict): Dictionary of sequences from the second dataset.

        Returns:
            Tuple[dict, List[str]]: A dictionary of common sequence mappings and
                                a list of comparison results.
        """
        common = {}

        for seq1, group1 in seq1_dict.items():
            if seq1 in seq2_dict:
                group2 = seq2_dict[seq1]
                # Check that the 2 groups have the same number of sequences
                if group1.count == group2.count:
                    if group1.count == 1:
                        common[group1.ids[0]] = group2.ids[0]
                    else:
                        self.comp.append(f"Matched 2 identical groups of sequences: {group1} and {group2}")
                        
                        # Map each ID in group1 to a possible group2 ID
                        possible_id2 = " OR ".join(group2.ids)
                        common.update({id1: possible_id2 for id1 in group1.ids})
                else:
                    self.comp.append(
                        f"Matched 2 different groups of sequences ({group1.count} vs {group2.count}): {group1} and {group2}"
                    )
        return common, self.comp
    
    def compare_seqs(self, fasta1 : Path, fasta2: Path) -> None:
        """
        Compare two FASTA files for common, unique, and differing sequences.

        Args:
            fasta1 : Path to the first FASTA file.
            fasta2 : Path to the second FASTA file.

        Result:
            A comparison log file to the output directory specified in `self.output_dir`.
        """

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
            self.comp.append(f"Sequences only in fasta1: {', '.join([str(ids) for ids in only1.values()])}")
        if only2:
            self.comp.append(f"Sequences only in fasta2: {', '.join([str(ids) for ids in only2.values()])}")

        if common:
            self.comp.append(f"Common ids: {', '.join([str(common_ids) for common_ids in common.items()])}")
        
        # Check for sequences with extra Ns
        self.compare_sequences_for_Ns(only1, only2)

        # Write results to file
        self.write_results()

    def write_results(self) -> None:
        """
        Write the comparison results to a file in the output directory.

        """
        output_file = Path.joinpath(self.output_dir, "compare.log")
        print(f"Writing results to {output_file}")
        with open(output_file, "w") as out_fh:
            for line in self.comp:
                out_fh.write(line + "\n")
    
    def compare_sequences_for_Ns(self, only1 : dict, only2: dict) -> None:
        """
        Compare sequences in `only1` and `only2` for differences in `N` content and length.

        Args:
            only1 (dict): Sequences unique to the first dataset, mapping sequence to group/identifier.
            only2 (dict): Sequences unique to the second dataset, mapping sequence to group/identifier.
        """
        # sequences which have extra N at the end
        if only1 and only2:
            for seq_1, name1 in only1.items():
                len1 = len(seq_1)
                seq1_N = seq_1.count("N")

                for seq_2, name2 in only2.items():
                    len2 = len(seq_2)
                    seq2_N = seq_2.count("N")

                    if abs(seq1_N - seq2_N) == abs(len1 - len2):
                        self.comp.append(f"Please check extra Ns added in your fasta2 in {name1} and {name2}")
                    else:
                        self.comp.append(
                            f"ALERT INSERTIONS at the end or diff assembly level {name1} and {name2}"
                        )

                    if len1 == len2:
                        if seq2_N > seq1_N:
                           self.comp.append(f"your fasta2 has more Ns, check {name1} and {name2}")
                        elif seq1_N > seq2_N:
                           self.comp.append(f"INSDC has more Ns, check {name1} and {name2}")
                        else:
                            self.comp.append(f"sequences have the same length, check {name1} and {name2}")

def parse_args(arg_list: list[str] | None) -> argparse.Namespace:
    """
    Parse command-line arguments for the genome sequence comparison tool.

    Args:
        arg_list (list[str] | None): A list of arguments to parse. If None, arguments
                                     are taken from sys.argv by default.

    Returns:
        argparse.Namespace: Parsed arguments as an argparse Namespace object.
    """
    parser = ArgumentParser(description="Compare sequences between two genomes")
    # Add filter arguments
    parser.add_argument("--fasta1_path", required=True, help="Path to INSDC fasta file")
    parser.add_argument("--fasta2_path", required=True, help="Path to user supplied fasta file")
    parser.add_argument("--output_dir", default=Path.cwd(), help="Directory to store the comparison report. Defaults to the current working directory.")
                        
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
    arg_list: Parsed arguments as an argparse Namespace object
    """
    args = parse_args(arg_list)

    CompareFasta(args.fasta1_path, args.fasta2_path, args.output_dir)

if __name__ == "__main__":
    main()