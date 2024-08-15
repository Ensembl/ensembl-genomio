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
    "check_chunk_size_and_tolerance",
    "chunk_fasta",
    "split_seq_by_chunk_size",
    "split_seq_by_n",
]


from contextlib import nullcontext
from io import TextIOWrapper
import logging
from os import PathLike
from pathlib import Path
import re
import sys
from typing import Callable, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def check_chunk_size_and_tolerance(
    chunk_size: int,
    chunk_tolerance: int,
    error_f: Callable[[str], None] = lambda x: sys.stderr.write(x) and sys.exit(-1) or sys.exit(-1),
):
    """Check the chunk size and the tolerance are positive
    and chunk size is not too small

    Args:
      chunk_size: Chunk size to check
      chunk_tolerance: Chunk tolerance to check

    Dies with` parser.error`
    """
    if chunk_size < 50_000:
        error_f(f"wrong '--chunk_size' value: '{chunk_size}'. should be greater then 50_000. exiting...")
    if chunk_tolerance < 0:
        error_f(f"wrong '--chunk_tolerance' value: '{chunk_tolerance}'. can't be less then 0. exiting...")


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


def chunk_fasta_stream(
    input_fasta: TextIOWrapper,
    chunk_size: int,
    chunk_size_tolerated: int,
    output_fasta: Optional[TextIOWrapper],
    individual_file_prefix: str,
    n_sequece_len: int = 0,
    chunk_sfx: str = "ens_chunk",
    append_offset_to_chunk_name: Optional[bool] = None,
    open_individual: Callable[[str], TextIOWrapper] = lambda name: open(name, "wt", encoding="utf-8"),
) -> list[str]:
    # output agp_lines list
    agp_lines = []

    # make sure not used for n_seq <= 0
    n_split_regex = None
    if n_sequece_len > 0:
        pattern = f"(N{{{n_sequece_len},}})"
        n_split_regex = re.compile(pattern)

    # process stream
    fasta_parser = SeqIO.parse(input_fasta, "fasta")
    for rec_count, rec in enumerate(fasta_parser, start=1):
        rec_name = str(rec.name)

        ends = split_seq_by_n(str(rec.seq), n_split_regex)
        ends = split_seq_by_chunk_size(ends, chunk_size, chunk_size_tolerated)

        offset = 0
        for chunk, chunk_end in enumerate(ends, start=1):
            chunk_name = f"{rec_name}_{chunk_sfx}_{chunk:03d}"
            chunk_file_name = f"{individual_file_prefix}.{rec_count:03d}.{chunk:03d}.fa"
            if append_offset_to_chunk_name:
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

            # if user specified chunks -- store each chunk in an individual output file
            with (
                output_fasta
                and nullcontext(output_fasta)  # no_type: ignore[arg-type]
                or open_individual(chunk_file_name)
            ) as out_file:
                out_file.write(tmp_record.format("fasta"))  # no_type: ignore[union-attr]

            del tmp_record
            offset = chunk_end

    return agp_lines


def get_tolerated_size(size: int, tolerance: int) -> int:
    if tolerance < 0:
        tolerance = 0

    tolerated_size = size
    tolerated_size += size * tolerance // 100
    return tolerated_size


def chunk_fasta(
    input_fasta_file: str,
    chunk_size: int,
    chunk_size_tolerated: int,
    out_file_name: str,
    individual_file_prefix: str,
    agp_output_file: Optional[str] = None,
    n_sequece_len: int = 0,
    chunk_sfx: str = "ens_chunk",
    append_offset_to_chunk_name: Optional[bool] = None,
):
    """ """

    # process input fasta
    with open_gz_file(input_fasta_file) as fasta:
        logging.info(
            f"splitting sequences from '{input_fasta_file}', chunk size {chunk_size:_}, \
                splitting on {n_sequece_len} Ns (0 -- disabled)"
        )
        # do not open a joined file if you plan to open many indvidual ones
        with (
            individual_file_prefix and nullcontext(None) or open(out_file_name, "wt", encoding="utf-8")
        ) as out_file_joined:
            agp_lines = chunk_fasta_stream(
                fasta,
                chunk_size,
                chunk_size_tolerated,
                out_file_joined,
                individual_file_prefix,
                n_sequece_len,
                chunk_sfx,
                append_offset_to_chunk_name,
            )

        # dump AGP
        if agp_output_file:
            with open(agp_output_file, "w", encoding="utf-8") as agp_out:
                agp_out.write("\n".join(agp_lines) + "\n")


def prepare_out_dir_for_individuals(dir_name: PathLike | str, file_part: str) -> Optional[str]:
    file_prefix = None
    if dir_name:
        dir_name.mkdir(parents=True, exist_ok=True)
        file_prefix = Path(dir_name, file_part)
    return file_prefix


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
        required=False,
        default=None,
        help="AGP file with chunks to contigs mapping.",
    )
    parser.add_argument(
        "--chunk_size",
        required=False,
        default=100_000_000,
        metavar="100_000_000",
        type=int,
        help="Maximum chunk size (should be greater then 50k).",
    )
    parser.add_argument(
        "--chunk_sfx",
        required=False,
        default="ens_chunk",
        type=str,
        help="Added to contig ID before chunk number.",
    )
    parser.add_argument(
        "--chunk_tolerance",
        required=False,
        default=0,
        type=int,
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
    parser.add_argument(
        "--add_offset",
        required=False,
        action="store_true",
        help="Append zero-based offset to chunk name ('_off_{offset}').",
    )

    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    check_chunk_size_and_tolerance(args.chunk_size, args.chunk_tolerance, error_f=parser.error)

    chunk_size_tolerated = get_tolerated_size(args.chunk_size, args.chunk_tolerance)

    file_prefix = prepare_out_dir_for_individuals(args.individual_out_dir, args.out or args.fasta_dna)

    chunk_fasta(
        args.fasta_dna,
        args.chunk_size,
        chunk_size_tolerated,
        args.out,
        individual_file_prefix=file_prefix,
        agp_output_file=args.agp_output_file,
        n_sequece_len=args.n_seq,
        chunk_sfx=args.chunk_sfx,
        append_offset_to_chunk_name=args.add_offset,
    )


if __name__ == "__main__":
    main()
