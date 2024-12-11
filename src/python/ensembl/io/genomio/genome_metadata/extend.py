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
"""Update a genome metadata file to include additional sequence regions (e.g. MT chromosome)."""

__all__ = [
    "get_additions",
    "get_gbff_regions",
    "get_report_regions_names",
    "amend_genome_metadata",
]

import csv
from os import PathLike
from pathlib import Path
import re
from typing import Dict, List, Tuple, Optional

from Bio import SeqIO

import ensembl.io.genomio
from ensembl.io.genomio.utils import get_json, print_json
from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


_VERSION_END = re.compile(r"\.\d+$")


def get_additions(report_path: PathLike, gbff_path: Optional[PathLike]) -> List[str]:
    """Returns all `seq_regions` that are mentioned in the report but that are not in the data.

    Args:
        report_path: Path to the report file.
        gbff_path: Path to the GBFF file.
    """
    gbff_regions = set(get_gbff_regions(gbff_path))
    report_regions = get_report_regions_names(report_path)
    additions = []
    for seq_region_name in report_regions:
        (genbank_seq_name, refseq_seq_name) = seq_region_name
        if genbank_seq_name not in gbff_regions and refseq_seq_name not in gbff_regions:
            if refseq_seq_name:
                additions.append(refseq_seq_name)
            else:
                additions.append(genbank_seq_name)
    additions = sorted(additions)
    return additions


def get_gbff_regions(gbff_path: Optional[PathLike]) -> List[str]:
    """Returns the `seq_region` data from a GBFF file.

    Args:
        gbff_path: GBFF file path to use.
    """
    seq_regions = []
    if gbff_path:
        with open_gz_file(gbff_path) as gbff_file:
            for record in SeqIO.parse(gbff_file, "genbank"):
                record_id = re.sub(_VERSION_END, "", record.id)
                seq_regions.append(record_id)
    return seq_regions


def _report_to_csv(report_path: PathLike) -> Tuple[str, Dict]:
    """Returns the assembly report as a CSV string, and its metadata as a dictionary.

    Args:
        report_path: Path to the assembly report file from INSDC/RefSeq.
    """
    data = ""
    metadata = {}
    with Path(report_path).open("r") as report:
        prev_line = ""
        for line in report:
            if line.startswith("#"):
                # Get metadata values if possible
                match = re.search(r"^#\s*([^:]+?):\s+(.+?)\s*$", line)
                if match:
                    metadata[match.group(1)] = match.group(2)
                prev_line = line
            else:
                if prev_line:
                    # Add previous line as header of CSV string, removing the initial "# "
                    data += prev_line[2:].strip() + "\n"
                    prev_line = ""
                data += line
    return data, metadata


def get_report_regions_names(report_path: PathLike) -> List[Tuple[str, str]]:
    """Returns a list of GenBank-RefSeq `seq_region` names from the assembly report file.

    Args:
        report_path: Path to the assembly report file from INSDC/RefSeq.
    """
    # Get the report in a CSV format, easier to manipulate
    report_csv, _ = _report_to_csv(report_path)
    # Feed the CSV string to the CSV reader
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


def amend_genome_metadata(
    genome_infile: PathLike,
    genome_outfile: PathLike,
    report_file: Optional[PathLike] = None,
    genbank_file: Optional[PathLike] = None,
) -> None:
    """
    Args:
        genome_infile: Genome metadata following the `src/python/ensembl/io/genomio/data/schemas/genome.json`.
        genome_outfile: Amended genome metadata file.
        report_file: INSDC/RefSeq sequences report file.
        genbank_file: INSDC/RefSeq GBFF file.
    """
    genome_metadata = get_json(genome_infile)
    # Get additional sequences in the assembly but not in the data
    if report_file:
        genbank_path = Path(genbank_file) if genbank_file else None
        additions = get_additions(report_file, genbank_path)
        if additions:
            genome_metadata["added_seq"] = {"region_name": additions}
    # Print out the file
    genome_outfile = Path(genome_outfile)
    print_json(genome_outfile, genome_metadata)


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--genome_infile",
        required=True,
        help="Input genome metadata file (following src/python/ensembl/io/genomio/data/schemas/genome.json)",
    )
    parser.add_argument_dst_path(
        "--genome_outfile", required=True, help="Path to the new amended genome metadata file"
    )
    parser.add_argument_src_path("--report_file", help="INSDC/RefSeq sequences report file")
    parser.add_argument_src_path("--genbank_file", help="INSDC/RefSeq GBFF file")
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    amend_genome_metadata(
        genome_infile=args.genome_infile,
        genome_outfile=args.genome_outfile,
        report_file=args.report_file,
        genbank_file=args.genbank_file,
    )
