#!env python3
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

import argparse
from pathlib import Path
from os import PathLike
from typing import List, Optional, Set

from Bio import SeqIO

from ensembl.io.genomio.utils.archive_utils import open_gz_file

exclude_seq_regions: List[str] = []


class GFFParserError(Exception):
    """Error while parsing a GFF file."""


def get_peptides_to_exclude(genbank_path: PathLike, seqr_to_exclude: Set[str]) -> Set[str]:
    """
    Extract peptide IDs from a genbank file that are in a given list of seq regions
    """
    peptides_to_exclude: Set[str] = set()
    with open_gz_file(genbank_path) as in_genbank:
        for record in SeqIO.parse(in_genbank, "genbank"):
            if record.id in seqr_to_exclude:
                print(f"Skip sequence {record.id}")
                for feat in record.features:
                    if feat.type == "CDS":
                        if "protein_id" in feat.qualifiers:
                            feat_id = feat.qualifiers["protein_id"]
                            peptides_to_exclude.add(feat_id[0])
                        else:
                            raise GFFParserError(f"Peptide without peptide ID ${feat}")
    return peptides_to_exclude


def prep_fasta_data(
    fasta_infile: PathLike,
    fasta_outfile: PathLike,
    peptide_mode: bool = False,
    genbank_infile: Optional[PathLike] = None,
) -> None:
    """
    Args:
        fasta_file: Input fasta file - DNA / Protein
        output_dir: Output folder for the fasta sequence file.
        peptide_mode: Process proteins not DNA
        genbank_infile: Input GenBank GBFF file
    """
    file_path = Path(fasta_infile)

    seqr_to_exclude = set(exclude_seq_regions)
    if peptide_mode:
        genbank_path = Path(genbank_infile)
        to_exclude = get_peptides_to_exclude(genbank_path, seqr_to_exclude)
    else:
        to_exclude = seqr_to_exclude

    # Copy and filter
    records = []

    # Final path
    try:
        with open_gz_file(file_path) as in_fasta:
            for record in SeqIO.parse(in_fasta, "fasta"):
                if record.id in to_exclude:
                    print(f"Skip record ${record.id}")
                else:
                    records.append(record)
    except FileNotFoundError:
        raise FileNotFoundError(f"{file_path} does not exist")

    with Path(fasta_outfile).open("w") as out_fasta:
        SeqIO.write(records, out_fasta, "fasta")


def main() -> None:
    """Module's entry-point."""
    parser = argparse.ArgumentParser(
        description=("Expand the genome metadata with information about the provider, assembly and gene"
                     "build version, and taxonomy."),
    )
    parser.add_argument("--fasta_infile", required=True, type=Path, help="Input FASTA file - DNA / Protein")
    parser.add_argument("--genbank_infile", type=Path, help="Input GenBank GBFF file (Optional)")
    parser.add_argument(
        "--fasta_outfile", type=Path, default=Path("output.fasta"), help="Output FASTA file path"
    )
    parser.add_argument(
        "--peptide_mode", type=bool, default=False, help="Process proteins not DNA (default: False)"
    )
    args = parser.parse_args()

    prep_fasta_data(**vars(args))
