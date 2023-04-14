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
"""TODO"""

from functools import partial
import gzip
from mimetypes import guess_type
from pathlib import Path
import argschema
from typing import Optional

from Bio import SeqIO

exclude_seq_regions = []


def peptides_to_exclude(genbank_path, seqr_to_exclude) -> dict():
    """
    Extract peptide IDs from a genbank file that are in a given list of seq regions
    """
    encoding = guess_type(genbank_path)[1]
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open
    peptides_to_exclude = dict()
    with _open(genbank_path) as in_genbank:
        for record in SeqIO.parse(in_genbank, "genbank"):
            if record.id in seqr_to_exclude:
                print("Skip sequence %s" % record.id)
                for feat in record.features:
                    if feat.type == "CDS":
                        if "protein_id" in feat.qualifiers:
                            feat_id = feat.qualifiers["protein_id"]
                            peptides_to_exclude[feat_id[0]] = True
                        else:
                            raise Exception(f"Peptide without peptide ID ${feat}")
    return peptides_to_exclude


def prep_fasta_data(
    fasta_infile: str,
    genebank_infile: str,
    peptide_mode: Optional[int] = 0,
) -> None:
    """
    Args:
        fasta_file: Input fasta file - DNA / Protein
        genbank_infile: Input genBank GBFF file (Optional)
        peptide_mode: Process proteins not DNA (Optional)
    """

    file_path = Path(fasta_infile)
    work_dir = Path("work_dir")
    seqr_to_exclude = exclude_seq_regions
    if peptide_mode:
        genbank_path = Path(genebank_infile)
        to_exclude = peptides_to_exclude(genbank_path, seqr_to_exclude)
    else:
        to_exclude = seqr_to_exclude

    work_dir.mkdir(parents=True, exist_ok=True)

    # Final file name
    new_file_name = fasta_infile + ".fa"
    # Final path
    final_path = work_dir / new_file_name
    # Copy and filter
    records = []
    encoding = guess_type(file_path)[1]
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open
    with _open(file_path) as in_fasta:
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
        required=True, metadata={"description": "Input fasta file - DNA / Protein"}
    )
    genbank_infile = argschema.fields.InputFile(
        required=False, metadata={"description": "Input genbank GBFF file (Optional)"}
    )
    peptide_mode = argschema.fields.Int(
        required=False, metadata={"description": "Process proteins not DNA (Optional)"}
    )


def main() -> None:
    """Module's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    prep_fasta_data(mod.args["fasta_infile"], mod.args["genbank_infile"], mod.args["peptide_mode"])
