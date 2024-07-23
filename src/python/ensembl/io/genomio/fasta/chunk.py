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
""""Split a set of nucleotide sequence(s) (.fasta, .gz) into smaller chunks."""

__all__ = [
    "split_seq_by_chunk_size",
    "split_seq_by_n",
]

import logging
from pathlib import Path
import re
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def split_seq_by_n(seq: str, split_pattern: re.Pattern) -> list:
    """Split a string into chunks at the positions where the
    pattern is found.

    The end point of each chunk will correspond to the end
    of the matching part.

    Args:
        seq: Sequence to be split into chunks.
        split_pattern: Pattern to search in the sequence.
    """
    seq_len = len(seq)
    if not split_pattern:
        return [seq_len]
    split_points = [m.end() for m in split_pattern.finditer(seq)]
    if not split_points or split_points[-1] != seq_len:
        split_points.append(seq_len)
    return split_points


def split_seq_by_chunk_size(ends: list[int], chunk_size: int, tolerated_size: Optional[int] = None) -> list:
    """Split list of end coordinates, to form chunks not longer then
    chunk_size.

    Args:
        ends: List of one or more chunk(s) to split a sequence.
        chunk_size: Size of chunks to split a sequence into.
        tolerated_size: Threshold to use instead of `chunk_size` to determine when to split a sequence.
    """
    if tolerated_size is None or tolerated_size < chunk_size:
        tolerated_size = chunk_size
    result = []
    offset = 0
    for chunk_end in ends:
        chunk_len = chunk_end - offset
        if chunk_len > tolerated_size:
            # exclude starts, as they are 0 or pushed as previous chunk_ends
            result += list(range(offset, chunk_end, chunk_size))[1:]
        result.append(chunk_end)
        offset = chunk_end
    return result


def main():
    """Module entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--fasta_dna",
        required=True,
        metavar="input[.fa | .gz]",
        help="Raw or compressed (.gz) FASTA file with DNA sequences to be split",
    )
    parser.add_argument(
        "--out",
        required=False,
        default="chunks.fna",
        help="Chunks output file",
    )
    parser.add_argument(
        "--individual_out_dir",
        required=False,
        default=None,
        type=Path,
        help="Output directory for writing files with individual chunks to. \
        If provided,`--out` value used as a filename prefix",
    )
    parser.add_argument_dst_path(
        "--agp_output_file",
        default=None,
        required=False,
        help="AGP file with chunks to contigs mapping.",
    )
    parser.add_argument(
        "--chunk_size",
        metavar="100_000_000",
        required=False,
        type=int,
        default=100_000_000,
        help="Maximum chunk size (should be greater then 50k).",
    )
    parser.add_argument(
        "--chunk_sfx",
        required=False,
        type=str,
        default="ens_chunk",
        help="Added to contig ID before chunk number.",
    )
    parser.add_argument(
        "--chunk_tolerance",
        required=False,
        type=int,
        default=0,
        help="Chunk size tolerance percentage. If the to-be-written chunk is longer \
            than the defined chunk size by less than the specified tolerance percentage,\
            it will not be split.",
    )
    parser.add_argument(
        "--n_seq",
        required=False,
        default=0,
        type=int,
        help="Split into chunks at positions of at least this number of N characters.",
    )
    parser.add_argument("--add_offset", required=False, action="store_true", help="Zero-based offset.")

    parser.add_log_arguments(add_log_file=True)

    args = parser.parse_args()

    if args.chunk_size < 50_000:
        parser.error(
            f"wrong '--chunk_size' value: '{args.chunk_size}'. should be greater then 50_000. exiting..."
        )
    if args.chunk_tolerance < 0:
        parser.error(
            f"wrong '--chunk_tolerance' value: '{args.chunk_tolerance}'. can't be less then 0. exiting..."
        )

    init_logging_with_args(args)

    tolerated_chunk_len = args.chunk_size
    tolerated_chunk_len += args.chunk_size * args.chunk_tolerance // 100

    # make sure not used for args.n_seq <= 0
    n_split_regex = None
    if args.n_seq > 0:
        pattern = f"(N{{{args.n_seq},}})"
        n_split_regex = re.compile(pattern)

    # output_file to write sequences to
    file_prefix = ""
    if args.individual_out_dir:
        if not args.out:
            args.out = args.fasta_dna
        args.individual_out_dir.mkdir(parents=True, exist_ok=True)
        file_prefix = Path(args.individual_out_dir, args.out)

    # process input fasta
    fasta_file = args.fasta_dna
    with open_gz_file(fasta_file) as fasta:
        agp_lines = []
        logging.info(
            f"splitting sequences from '{fasta_file}', chunk size {args.chunk_size:_}, \
                splitting on {args.n_seq} Ns (0 -- disabled)"
        )

        fasta_parser = SeqIO.parse(fasta, "fasta")
        for rec_count, rec in enumerate(fasta_parser, start=1):
            rec_name = str(rec.name)

            ends = split_seq_by_n(str(rec.seq), n_split_regex)
            ends = split_seq_by_chunk_size(ends, args.chunk_size, tolerated_chunk_len)

            offset = 0
            for chunk, chunk_end in enumerate(ends, start=1):
                chunk_name = f"{rec_name}_{args.chunk_sfx}_{chunk:03d}"
                chunk_file_name = f"{file_prefix}.{rec_count:03d}.{chunk:03d}.fa"
                if args.add_offset:
                    chunk_name += f"_off_{offset}"

                rec_from = offset + 1
                rec_to = chunk_end
                chunk_len = chunk_end - offset

                # form agp lines
                agp_line = f"{rec_name}\t{rec_from}\t{rec_to}\t{chunk}\tW\t{chunk_name}\t1\t{chunk_len}\t+"
                agp_lines.append(agp_line)

                # use agp lines as fasta description
                agp_line = agp_line.replace("\t", " ")
                logging.info(f"Dumping {chunk_name} AGP {agp_line}")

                # get slice and put it out
                tmp_record = SeqRecord(
                    Seq(rec.seq[offset:chunk_end]), id=chunk_name, description=f"AGP {agp_line}", name=""
                )

                # if user specified chunks stored in individual output files
                if args.individual_out_dir:
                    with open(chunk_file_name, "wt", encoding="utf-8") as out_file:
                        out_file.write(tmp_record.format("fasta"))
                else:
                    with open(args.out, "wt", encoding="utf-8") as out_file:
                        out_file.write(tmp_record.format("fasta"))

                del tmp_record
                offset = chunk_end

        # dump AGP
        if args.agp_output_file:
            with open(args.agp_output_file, "w", encoding="utf-8") as agp_out:
                agp_out.write("\n".join(agp_lines))

        if not args.individual_out_dir:
            out_file.close()
