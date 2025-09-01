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
"""Compares the INSDC and core fasta files based on MD5sum and IDs."""

import argparse
import logging
from pathlib import Path
import re
from typing import Any, Dict, List, Tuple

from Bio import SeqIO

import ensembl.io.genomio
from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

"""Compare sequences between two genomes"""

__all__ = [
    "SeqGroup",
    "CompareFasta",
]


class SeqGroup:
    """Represents a group of sequence identifiers and maintains a count of them."""

    def __init__(self, identifier: str | None = None) -> None:
        """
        Initializes a SeqGroup instance.

        Args:
            identifier (str | None, optional): The first identifier to add to the group.
                                               Defaults to None, which adds "None" as the identifier.
        """
        self.ids: List[str] = []
        if identifier:
            self.add_id(identifier)
        else:
            self.add_id("None")
        self.count = len(self.ids)

    def __str__(self) -> str:
        """
        Returns a comma-separated string of sequence identifiers.

        Returns:
            str: A string representation of the sequence group.
        """
        return ", ".join(self.ids)

    def add_id(self, identifier: str | None = None) -> None:
        """
        Adds a new identifier to the group and updates the count.

        Args:
            identifier (str | None, optional): The identifier to add. If None, "None" is added instead.
        """
        self.ids.append(identifier if identifier else "None")
        self.count = len(self.ids)


class CompareFasta:
    """Read and compare the fasta sequences"""

    def __init__(self, fasta1: Path, fasta2: Path, output_dir: str) -> None:
        self.fasta1 = Path(fasta1)
        self.fasta2 = Path(fasta2)
        self.output_dir = Path(output_dir)
        self.comp: List[str] = []

    def compare_seqs(self) -> None:
        """
        Compare two FASTA files for common, unique, and differing sequences.
        """
        seq1 = self.read_fasta(self.fasta1)
        seq2 = self.read_fasta(self.fasta2)

        # Compare sequences
        seq1_dict = self.build_seq_dict(seq1)
        seq2_dict = self.build_seq_dict(seq2)

        # Compare number of sequences
        if len(seq1) != len(seq2):
            self.comp.append(
                "WARNING: Different number of sequences: "
                f"fasta1 [ n = {len(seq1)} ] -Vs- fasta2 [ n = {len(seq2)} ]"
            )
            logging.warning("Different number of sequences: fasta1 compared to fasta2")

        common, group_comp = self.find_common_groups(seq1_dict, seq2_dict)
        self.comp += group_comp

        # Sequences that are not common
        only1 = {seq: group for seq, group in seq1_dict.items() if not seq in seq2_dict}

        only2 = {seq: group for seq, group in seq2_dict.items() if not seq in seq1_dict}

        if only1:
            self.comp.append(f"Sequences only in Fasta_1: {', '.join([str(ids) for ids in only1.values()])}")
        if only2:
            self.comp.append(f"Sequences only in Fasta_2: {', '.join([str(ids) for ids in only2.values()])}")

        if common:
            self.comp.append(f"Common ids: {', '.join([str(common_ids) for common_ids in common.items()])}")

        # Check for sequences with extra Ns
        if only1 and only2:
            self.compare_seq_for_Ns(only1, only2)

        # Write results to file
        self.write_results()

    def read_fasta(self, fasta_path: Path) -> dict:
        """
        Reads a FASTA file and returns a dictionary mapping sequence IDs to sequences.

        Args:
            fasta_path (Path): Path to the FASTA file. Supports gzipped files.

        Returns:
            dict: A dictionary where keys are sequence IDs and values are sequences
              with all non-CGTA characters replaced by 'N'.
        """
        logging.info(f"Read fasta file {fasta_path}")
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
        seqs_dict: Dict[str, SeqGroup] = {}
        for name, seq in seqs.items():
            if seq in seqs_dict:
                seqs_dict[seq].add_id(name)
            else:
                seqs_dict[seq] = SeqGroup(name)

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
                        "Matched 2 different groups of sequences"
                        f" ({group1.count} vs {group2.count}): {group1} and {group2}"
                    )

        return common, self.comp

    def write_results(self) -> None:
        """
        Write the comparison results to a file in the output directory.

        """
        output_file = Path.joinpath(self.output_dir, "compare.log")
        observed_compare = set()

        logging.info(f"Writing results to {output_file}")
        with open(output_file, "w") as out_fh:
            for line in self.comp:
                if line not in observed_compare:
                    observed_compare.add(line)
                    out_fh.write(str(line) + "\n")

    def compare_seq_for_Ns(self, only1: dict, only2: dict) -> None:
        """
        Compare sequences in `only1` and `only2` for differences in `N` content and length.

        Args:
            only1 (dict): Sequences unique to the first dataset, mapping sequence to group/identifier.
            only2 (dict): Sequences unique to the second dataset, mapping sequence to group/identifier.
        """
        # sequences which have extra N at the end
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
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_argument("--fasta1", required=True, help="Path to INSDC fasta file")
    parser.add_argument("--fasta2", required=True, help="Query  fasta to check against INSDC fasta.")
    parser.add_argument(
        "--output_dir",
        default=Path.cwd(),
        help="Directory to store the comparison report. Defaults to the current working directory.",
    )
    # Add flags
    parser.add_argument(
        "--compare_seq_region",
        action="store_true",
        help="Enable compare seq_region mode, i.e. use seq_region for sequence comparison",
    )
    parser.add_log_arguments()
    return parser.parse_args(arg_list)


def main(arg_list: list[str] | None = None) -> None:
    """Main script entry-point.

    Args:
    arg_list: Parsed arguments as an argparse Namespace object
    """
    args = parse_args(arg_list)
    init_logging_with_args(args)
    compare_initialise = CompareFasta(args.fasta1, args.fasta2, args.output_dir)
    # Perform the comparison explicitly
    compare_initialise.compare_seqs()
