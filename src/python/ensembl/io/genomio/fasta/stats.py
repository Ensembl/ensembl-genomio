# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Compute longest, total, and count statistics for FASTA records."""

import argparse
from dataclasses import dataclass
import logging
from pathlib import Path

import ensembl.io.genomio
from ensembl.io.genomio.utils import json_utils
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.archive import open_gz_file
from ensembl.utils.logging import init_logging_with_args

__all__ = [
    "FastaStats",
    "compute_fasta_stats",
]


@dataclass(frozen=True)
class FastaStats:
    """Basic FASTA sequence length summary."""

    longest: int
    total: int
    n_seqs: int


def compute_fasta_stats(fasta_file: Path, output_file: Path | None) -> None:
    """
    Compute basic FASTA stats in a single streaming pass.

    Args:
        fasta_file: Path to raw or compressed input FASTA file.
        output_file: Path to the output stats JSON file (defaults to ``<fasta_file>.stats.json`` if `None`).
    """
    longest = 0
    total = 0
    n_seqs = 0
    current = 0

    with open_gz_file(fasta_file) as fh:
        for raw_line in fh:
            if raw_line.startswith(">"):
                if n_seqs > 0:
                    longest = max(longest, current)
                    total += current
                n_seqs += 1
                current = 0
                continue

            line = raw_line.strip()
            current += len(line)

    if n_seqs > 0:
        longest = max(longest, current)
        total += current

    stats = FastaStats(longest=longest, total=total, n_seqs=n_seqs)
    json_content = {
        "longest_seq": stats.longest,
        "total_seq_length": stats.total,
        "nr_seqs": stats.n_seqs,
    }

    output_file = output_file or Path(fasta_file).with_suffix(".stats.json")
    json_utils.print_json(output_file, json_content)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Args:
        argv: Optional command-line arguments. If `None`, argparse reads from ``sys.argv``.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--fasta",
        required=True,
        help="Input FASTA, optionally gzipped.",
    )
    parser.add_argument_dst_path(
        "--output",
        metavar="STATS.json",
        default=argparse.SUPPRESS,
        help="Output stats JSON file.",
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()

    args = parser.parse_args(argv)
    init_logging_with_args(args)

    return args


def main(argv: list[str] | None = None) -> None:
    """
    Run the FASTA stats command-line interface.

    Args:
        argv: Optional command-line arguments. If `None`, argparse reads from ``sys.argv``.
    """
    args = parse_args(argv)

    try:
        compute_fasta_stats(
            fasta_file=args.fasta,
            output_file=getattr(args, "output", None),
        )
    except Exception:
        logging.exception(f"Error processing FASTA file {args.fasta}")
        raise
