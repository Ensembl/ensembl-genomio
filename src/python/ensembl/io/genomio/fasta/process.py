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
"""Takes a FASTA file (DNA or peptide), cleans it up and optionally excludes some IDs."""

__all__ = ["FastaParserError", "get_peptides_to_exclude", "prep_fasta_data"]

import logging
from pathlib import Path
from os import PathLike
from typing import List, Optional, Set

from Bio import SeqIO

import ensembl.io.genomio
from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


exclude_seq_regions: List[str] = []


class FastaParserError(Exception):
    """Error while parsing a FASTA file."""


def get_peptides_to_exclude(genbank_path: PathLike, seqr_to_exclude: Set[str]) -> Set[str]:
    """
    Extract peptide IDs from a genbank file that are in a given list of seq regions
    """
    peptides_to_exclude: Set[str] = set()
    with open_gz_file(genbank_path) as in_genbank:
        for record in SeqIO.parse(in_genbank, "genbank"):
            if record.id in seqr_to_exclude:
                logging.info(f"Skip sequence {record.id}")
                for feat in record.features:
                    if feat.type == "CDS":
                        if "protein_id" in feat.qualifiers:
                            feat_id = feat.qualifiers["protein_id"]
                            peptides_to_exclude.add(feat_id[0])
                        else:
                            raise FastaParserError(f"Peptide without peptide ID ${feat}")
    return peptides_to_exclude


def prep_fasta_data(
    fasta_infile: PathLike,
    genbank_infile: Optional[PathLike],
    fasta_outfile: PathLike,
    peptide_mode: bool = False,
) -> None:
    """
    Args:
        fasta_file: Input FASTA file - DNA / Protein
        genbank_infile: Input GenBank GBFF file (Optional)
        fasta_outfile: Output FASTA sequence file.
        peptide_mode: Process proteins instead of DNA
    """
    file_path = Path(fasta_infile)

    to_exclude = set()
    seqr_to_exclude = set(exclude_seq_regions)
    if peptide_mode:
        if genbank_infile is not None:
            genbank_path = Path(genbank_infile)
            to_exclude = get_peptides_to_exclude(genbank_path, seqr_to_exclude)
    else:
        to_exclude = seqr_to_exclude

    # Copy and filter
    records = []

    # Final path
    with open_gz_file(file_path) as in_fasta:
        for record in SeqIO.parse(in_fasta, "fasta"):
            if record.id in to_exclude:
                logging.info(f"Skip record ${record.id}")
            else:
                records.append(record)
    with Path(fasta_outfile).open("w") as out_fasta:
        SeqIO.write(records, out_fasta, "fasta")


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description="Clean-up a given FASTA file to remove unwanted elements.")
    parser.add_argument_src_path("--fasta_infile", required=True, help="Input FASTA file - DNA / Protein")
    parser.add_argument_src_path("--genbank_infile", help="Input GenBank GBFF file")
    parser.add_argument_dst_path("--fasta_outfile", required=True, help="Output FASTA file")
    parser.add_argument("--peptide_mode", action="store_true", help="Process proteins instead of DNA")
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    prep_fasta_data(
        fasta_infile=args.fasta_infile,
        genbank_infile=args.genbank_infile,
        fasta_outfile=args.fasta_outfile,
        peptide_mode=args.peptide_mode,
    )
