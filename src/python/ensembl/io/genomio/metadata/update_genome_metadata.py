#!/usr/bin/env python
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
"""Add more metadata to the genome metadata file, including added seq_regions (e.g. MT chromosome)."""

import csv
from os import PathLike
from pathlib import Path
import re
from typing import List, Tuple, Optional
import argschema

from Bio import SeqIO

from ensembl.io.genomio.utils import print_json, get_json
from ensembl.io.genomio.utils.archive_utils import open_gz_file


class MissingDataError(Exception):
    """Used if some data is missing from the report file."""

    def __init__(self, report_path: PathLike, accession: str, msg: str):
        report_msg = f"Can't get data for {accession} in report {report_path}"
        if msg:
            report_msg = f"{report_msg}: {msg}"
        self.msg = report_msg


def get_additions(report_path: Path, gbff_path: Path) -> List[str]:
    """Returns all `seq_regions` that are mentioned in the report but that are not in the data.

    Args:
        report_path: Path to the report file.
        gbff_path: Path to the GBFF file.
    """
    gbff_regions = set(get_gbff_regions(gbff_path))
    report_regions = set(get_report_regions_names(report_path))
    additions = list(report_regions.difference(gbff_regions))
    return additions


def get_gbff_regions(gbff_path: Path) -> List[str]:
    """Returns the `seq_region` data from the GBFF file.

    Args:
        gbff_path: Gbff file path to use.
    """
    if not gbff_path:
        return []

    seq_regions = []
    with open_gz_file(gbff_path) as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            seq_regions.append(record.id)
    return seq_regions


def _report_to_csv(report_path: Path) -> Tuple[str, dict]:
    """Returns an assembly report as a CSV string, and the head metadata as a dict.

    Args:
        report_path: Path to a `seq_region` file from INSDC/RefSeq.

    """
    with report_path.open("r") as report:
        data = ""
        metadata = {}
        last_head = ""
        for line in report:
            # Ignore header
            if line.startswith("#"):
                # Get metadata values if possible
                match = re.search("# (.+?): (.+?)$", line)
                if match:
                    metadata[match.group(1)] = match.group(2)
                last_head = line
                continue
            if last_head:
                data += last_head[2:].strip() + "\n"
                last_head = ""
            data += line
        return data, metadata


def get_report_regions_names(report_path: Path) -> List[str]:
    """Returns a list of `seq_region` names from the report file.

    Args:
        report_path: Path to the seq_regions report from INSDC/RefSeq.
    """
    # Get the report in a CSV format, easier to manipulate
    report_csv, _ = _report_to_csv(report_path)

    # Feed the csv string to the CSV reader
    reader = csv.DictReader(report_csv.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE)

    # Create the seq_regions
    seq_regions = []
    for row in reader:
        refseq_name = row["RefSeq-Accn"]
        genbank_name = row["GenBank-Accn"]
        name = ""
        if genbank_name and genbank_name != "na":
            name = genbank_name
        elif refseq_name and refseq_name != "na":
            name = refseq_name
        if name:
            seq_regions.append(name)
    return seq_regions


def amend_genomic_metadata(
    genome_infile: PathLike,
    insdc_refseq_report_infile: PathLike,
    genbank_infile: PathLike,
    output_dir: PathLike,
    brc4_mode: Optional[int] = 1,
) -> Path:
    """
    Args:
        genome_infile: Genome data following the schemas/genome_schema.json.
        INSDC_RefSeq_report_infile: Path to the INSDC/RefSeq sequences report to parse.
        genbank_infile: Path to the INSDC/RefSeq gbff file to parse.
        output_dir: Directory where the new file will be created.
        brc4_mode: Activate BRC4 mode (default).
    """

    # Create dedicated work dir
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Make dict from Genome JSON
    genome_metadata = get_json(genome_infile)

    # Final file name
    metadata_type = "genome"
    new_file_name = f"{metadata_type}_amended.json"
    final_path = output_dir / new_file_name
    # use_refseq = self.param("accession").startswith("GCF_")

    # Get additional sequences in the assembly but not in the data
    additions = get_additions(Path(insdc_refseq_report_infile), Path(genbank_infile))
    if additions:
        genome_metadata["added_seq"] = {"region_name": additions}

    # Possible brc specific
    if brc4_mode:
        pass

    # Print out the file
    print_json(final_path, genome_metadata)

    return final_path


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by the entry point of this module."""

    genome_infile = argschema.fields.InputFile(
        required=True, metadata={"description": "Input Genome data. Following the schemas/genome_schema.json"}
    )
    INSDC_RefSeq_report_infile = argschema.fields.InputFile(
        required=False, metadata={"description": "Path to the INSDC/RefSeq sequences report to parse"}
    )
    genbank_infile = argschema.fields.InputFile(
        required=False, metadata={"description": "Path to the INSDC/RefSeq gbff file to parse"}
    )
    output_dir = argschema.fields.OutputDir(
        required=False,
        dump_default=".",
        metadata={"description": "Directory where the new amended genome file will be created"},
    )
    brc4_mode = argschema.fields.Int(required=False, metadata={"description": "Activate BRC4 mode (default)"})


def main() -> None:
    """Module's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    amended_genome_file = amend_genomic_metadata(
        mod.args["genome_infile"],
        mod.args["INSDC_RefSeq_report_infile"],
        mod.args["genbank_infile"],
        mod.args["output_dir"],
        mod.args["brc4_mode"],
    )
    # Flow out the file and type
    output = {"metadata_type": "genome", "metadata_json": str(amended_genome_file)}
    print_json(mod.args["output_json"], output)
