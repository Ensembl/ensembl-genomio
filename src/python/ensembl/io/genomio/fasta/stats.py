#!/usr/bin/env python3

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
import logging
from dataclasses import dataclass
from pathlib import Path

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.archive import open_gz_file
from ensembl.utils.logging import init_logging_with_args

__all__ = [
    "FastaStats",
    "compute_fasta_stats",
    "parse_args",
    "main",
]


@dataclass(frozen=True)
class FastaStats:
    longest: int
    total: int
    n_seqs: int


def _write_fasta_stats(stats: FastaStats, output_file: Path) -> None:
    """
    Write FASTA statistics to the output file.

    Args:
        stats: FASTA statistics to write.
        output: Path to the output text file.
    """
    output_file.write_text(
        f"{stats.longest} {stats.total} {stats.n_seqs}\n",
        encoding="utf-8",
    )


def compute_fasta_stats(fasta_file: Path, output_file: Path | None) -> None:
    """
    Compute basic FASTA stats in a single streaming pass.

    Args:
        fasta_file: Path to raw or compressed input FASTA file.
        output_file: Path to the output stats file.

    Returns:
        None
    """
    longest = 0
    total = 0
    n_seqs = 0
    current = 0
    saw_record = False

    with open_gz_file(fasta_file) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if line.startswith(">"):
                if saw_record:
                    longest = max(longest, current)
                    total += current
                n_seqs += 1
                current = 0
                saw_record = True
                continue

            current += len(line.strip())

    if saw_record:
        longest = max(longest, current)
        total += current

    output_file = output_file or Path(fasta_file).with_suffix(".stats.txt")
    _write_fasta_stats(FastaStats(longest=longest, total=total, n_seqs=n_seqs), output_file)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Args:
        argv: Optional command-line arguments. If None, argparse reads from sys.argv.

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
        metavar="STATS.txt",
        default=argparse.SUPPRESS,
        help="Output stats text file.",
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
        argv: Optional command-line arguments. If None, argparse reads from sys.argv.

    Returns:
        int: Exit status code.
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
