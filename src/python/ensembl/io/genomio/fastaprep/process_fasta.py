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
"""Takes a fasta file (DNA or peptide), clean it up and optionally exclude some ids."""

from mimetypes import guess_type
from pathlib import Path
from os import PathLike
from typing import List, Optional, Set

import argschema

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
    genebank_infile: PathLike,
    output_dir: PathLike,
    peptide_mode: Optional[int] = 0,
) -> None:
    """
    Args:
        fasta_file: Input fasta file - DNA / Protein
        genbank_infile: Input genBank GBFF file (Optional)
        output_dir: Output folder for the fasta sequence file.
        peptide_mode: Process proteins not DNA (Optional)
    """
    file_path = Path(fasta_infile)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    seqr_to_exclude = set(exclude_seq_regions)
    if peptide_mode:
        genbank_path = Path(genebank_infile)
        to_exclude = get_peptides_to_exclude(genbank_path, seqr_to_exclude)
    else:
        to_exclude = seqr_to_exclude

    # Final file name
    new_file_name = file_path.stem

    # Copy and filter
    records = []
    encoding = guess_type(file_path)[1]

    if encoding == "gzip":
        new_file_name = Path(new_file_name).stem

    # Final path
    final_path = output_dir / f"{new_file_name}.fa"
    with open_gz_file(file_path) as in_fasta:
        for record in SeqIO.parse(in_fasta, "fasta"):
            if record.id in to_exclude:
                print(f"Skip record ${record.id}")
            else:
                records.append(record)
    with final_path.open("w") as out_fasta:
        SeqIO.write(records, out_fasta, "fasta")


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by the entry point of this module."""

    fasta_infile = argschema.fields.InputFile(
        required=True,
        metadata={"description": "Input fasta file - DNA / Protein"},
    )
    genbank_infile = argschema.fields.InputFile(
        required=False,
        metadata={"description": "Input genbank GBFF file (Optional)"},
    )
    output_dir = argschema.fields.OutputDir(
        required=False,
        dump_default=".",
        metadata={"description": "Dir name where processed fasta data will reside"},
    )
    peptide_mode = argschema.fields.Int(
        required=False,
        metadata={"description": "Process proteins not DNA (Optional)"},
    )


def main() -> None:
    """Module's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    prep_fasta_data(
        mod.args["fasta_infile"], mod.args["genbank_infile"], mod.args["output_dir"], mod.args["peptide_mode"]
    )
