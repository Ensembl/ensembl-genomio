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

__all__ = [
    "MissingDataError",
    "get_additions",
    "get_gbff_regions",
    "get_report_regions_names",
    "amend_genomic_metadata",
]

import csv
from os import PathLike
from pathlib import Path
import re
from typing import List, Tuple, Optional

from Bio import SeqIO

from ensembl.io.genomio.utils import get_json, open_gz_file, print_json
from ensembl.utils.argparse import ArgumentParser


_VERSION_END = re.compile(r"\.\d+$")


class MissingDataError(Exception):
    """Used if some data is missing from the report file."""

    def __init__(self, report_path: PathLike, accession: str, msg: str):
        report_msg = f"Can't get data for {accession} in report {report_path}"
        if msg:
            report_msg = f"{report_msg}: {msg}"
        self.msg = report_msg


def get_additions(report_path: Path, gbff_path: Optional[Path]) -> List[str]:
    """Returns all `seq_regions` that are mentioned in the report but that are not in the data.

    Args:
        report_path: Path to the report file.
        gbff_path: Path to the GBFF file.
    """
    gbff_regions = set(get_gbff_regions(gbff_path))
    report_regions = get_report_regions_names(report_path)

    additions = []
    for rep_seq in report_regions:
        (rs_seq, gb_seq) = rep_seq
        if rs_seq not in gbff_regions and gb_seq not in gbff_regions:
            if rs_seq:
                additions.append(rs_seq)
            else:
                additions.append(gb_seq)
    additions = sorted(additions)
    return additions


def get_gbff_regions(gbff_path: Optional[Path]) -> List[str]:
    """Returns the `seq_region` data from the GBFF file.

    Args:
        gbff_path: Gbff file path to use.
    """
    if not gbff_path:
        return []

    seq_regions = []
    with open_gz_file(gbff_path) as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            record_id = re.sub(_VERSION_END, "", record.id)
            seq_regions.append(record_id)
    return seq_regions


def _report_to_csv(report_path: Path) -> Tuple[str, dict]:
    """Returns an assembly report as a CSV string, and the head metadata as a dict.

    Args:
        report_path: Path to a `seq_region` file from INSDC/RefSeq.

    """
    data = ""
    metadata = {}
    with report_path.open("r") as report:
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


def get_report_regions_names(report_path: Path) -> List[Tuple[str, str]]:
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

        if refseq_name == "na":
            refseq_name = ""
        if genbank_name == "na":
            genbank_name = ""
        refseq_name = re.sub(_VERSION_END, "", refseq_name)
        genbank_name = re.sub(_VERSION_END, "", genbank_name)
        seq_regions.append((genbank_name, refseq_name))
    return seq_regions


def amend_genomic_metadata(
    genome_infile: PathLike,
    genome_outfile: PathLike,
    report_file: Optional[PathLike] = None,
    genbank_infile: Optional[PathLike] = None,
) -> None:
    """
    Args:
        genome_infile: Genome data following the schemas/genome_schema.json.
        genome_outfile: Amended genome data file.
        report_file: INSDC/RefSeq sequences report file.
        genbank_infile: INSDC/RefSeq GBFF file.
    """
    genome_metadata = get_json(genome_infile)

    # Get additional sequences in the assembly but not in the data
    if report_file:
        gbff_path = Path(genbank_infile) if genbank_infile else None
        additions = get_additions(Path(report_file), gbff_path)
        if additions:
            genome_metadata["added_seq"] = {"region_name": additions}

    # Print out the file
    genome_outfile = Path(genome_outfile)
    print_json(genome_outfile, genome_metadata)


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(
        description="Update genome metadata file to include additional sequence regions (e.g. MT chromosome)."
    )
    parser.add_argument_src_path(
        "--genome_infile", required=True, help="Input genome file (following the schemas/genome_schema.json)"
    )
    parser.add_argument_dst_path(
        "--genome_outfile", required=True, help="Path to the new amended genome metadata file"
    )
    parser.add_argument_src_path("--report_file", help="INSDC/RefSeq sequences report file")
    parser.add_argument_src_path("--genbank_infile", help="INSDC/RefSeq GBFF file")
    args = parser.parse_args()

    amend_genomic_metadata(**vars(args))
