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
"""Split a set of nucleotide sequence(s) (.fasta, .gz) into smaller chunks."""

__all__ = [
    "check_chunk_size_and_tolerance",
    "chunk_fasta",
    "chunk_fasta_stream",
    "get_tolerated_size",
    "prepare_out_dir_for_individuals",
    "split_seq_by_chunk_size",
    "split_seq_by_n",
]


from contextlib import nullcontext
from io import TextIOWrapper
import logging
from pathlib import Path
import re
from typing import Any, Callable, ContextManager, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import ensembl.io.genomio
from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def _on_value_error(msg: str) -> None:
    """A default on_value_handler, raises ValueError with the provided `msg`.

    Args:
        msg: A message to raise ValueError with.
    """
    raise ValueError(msg)


def check_chunk_size_and_tolerance(
    chunk_size: int,
    chunk_tolerance: int,
    error_f: Callable[[str], None] = _on_value_error,
) -> None:
    """Check the chunk size and the tolerance are positive and chunk size is not too small

    Args:
        chunk_size: Chunk size to check
        chunk_tolerance: Chunk tolerance to check

    Dies:
        If checks failed dies with `parser.error`
    """
    if chunk_size < 50_000:
        error_f(f"wrong '--chunk_size' value: '{chunk_size}'. should be greater then 50_000. exiting...")
    if chunk_tolerance < 0:
        error_f(f"wrong '--chunk_tolerance' value: '{chunk_tolerance}'. can't be less then 0. exiting...")


def split_seq_by_n(seq: str, split_pattern: Optional[re.Pattern]) -> list[int]:
    """Split a string into chunks at the positions where the
    pattern is found. `N`s (pattern) are appended to the chunk on the left.

    The end point of each chunk will correspond to the end
    of the matching part.

    Args:
        seq: Sequence to be split into chunks.
        split_pattern: Pattern to search in the sequence.

    Returns:
        List with open coordinates of the chunks ends (or with only a single sequence length).
    """
    seq_len = len(seq)
    if not split_pattern:
        return [seq_len]
    split_points = [m.end() for m in split_pattern.finditer(seq)]
    if not split_points or split_points[-1] != seq_len:
        split_points.append(seq_len)
    return split_points


def split_seq_by_chunk_size(
    ends: list[int], chunk_size: int, tolerated_size: Optional[int] = None
) -> list[int]:
    """Split list of end coordinates, to form chunks not longer then
    chunk_size.

    Args:
        ends: List of one or more chunk(s) to split a sequence.
        chunk_size: Size of chunks to split a sequence into.
        tolerated_size: Threshold to use instead of `chunk_size` to determine when to split a sequence.

    Returns:
        List with open coordinates of the chunks ends (or with only a single sequence length).
    """
    if not ends or chunk_size < 1:
        return ends

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


def _individual_file_opener(name: str) -> TextIOWrapper:
    """Is used to open file for an individual chunk sequence.

    Args:
        name: Name of the file to open
    """
    return open(name, "wt", encoding="utf-8")


def chunk_fasta_stream(
    input_fasta: TextIOWrapper,
    chunk_size: int,
    chunk_size_tolerated: int,
    output_fasta: Optional[TextIOWrapper] | nullcontext[Any],
    individual_file_prefix: Optional[Path],
    *,
    n_sequence_len: int = 0,
    chunk_sfx: str = "ens_chunk",
    append_offset_to_chunk_name: Optional[bool] = None,
    open_individual: Callable[[str], ContextManager[Any]] = _individual_file_opener,
) -> list[str]:
    """Split input TextIOWrapper stream with fasta into a smaller chunks based on
    stretches of "N"s and then based on chunk_size_tolerated and store either to
    the output_fasta stream (if valid) or to the files created by
    invocation of the `open_individual` callable.

    Args:
        input_fasta: Input FASTA as the TextIOWrapper stream.
        chunk_size: Size of the chunks to split into.
        chunk_size_tolerated: If more flexibility allowed, use this as the maximum size of a chunk.
        output_fasta: Output FASTA TextIOWrapper stream to store the chunks into,
                if none or False, `open_individual` is used (see below).
        individual_file_prefix: A file path prefix including dirs and filenames part to use as a
                first part of the chunk file name.
        n_sequence_len: Length of the stretch of `N`s to split at; not slitting on `N`s if 0.
        chunk_sfx: A string to put between the original sequence name and the chunk suffix.
        append_offset_to_chunk_name: A flag to append 0-based offset (`_off_{offset}`) to the chunk name.
        open_individual: A callable taking filename as an input to generate the output file for
                individual contig if out_put FASTA is `false` of `None`, folders should be preexisting.
    """

    chunk_size_tolerated = max(chunk_size, chunk_size_tolerated)
    # output agp_lines list
    agp_lines = []

    # make sure not used for n_seq <= 0
    n_split_regex = None
    if n_sequence_len > 0:
        pattern = f"(N{{{n_sequence_len},}})"
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
            chunk_file_name = ""
            if individual_file_prefix:
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
                Seq(rec.seq[offset:chunk_end]),
                id=chunk_name,
                description=f"AGP {agp_line}",
                name="",
            )

            # if user specified chunks -- store each chunk in an individual output file
            with output_fasta and nullcontext(output_fasta) or open_individual(chunk_file_name) as out_file:
                out_file.write(tmp_record.format("fasta"))  # type: ignore[union-attr]

            del tmp_record
            offset = chunk_end

    return agp_lines


def get_tolerated_size(size: int, tolerance: int) -> int:
    """Calculate max tolerated size of the chunk based on initial size and percent of allowed deviation.

    Args:
        size: Base chunk size
        tolerance: Percent of allowed deviance as integer.

    Returns:
        Maximum tolerated chunk size.
    """
    tolerance = max(tolerance, 0)

    tolerated_size = size
    tolerated_size += size * tolerance // 100
    return tolerated_size


def chunk_fasta(
    input_fasta_file: str,
    chunk_size: int,
    chunk_size_tolerated: int,
    out_file_name: str,
    individual_file_prefix: Optional[Path],
    *,
    agp_output_file: Optional[str] = None,
    n_sequence_len: int = 0,
    chunk_sfx: str = "ens_chunk",
    append_offset_to_chunk_name: Optional[bool] = None,
) -> None:
    """Open `input_fasta_file` and split into a smaller chunks based on
    stretches of "N"s and then based on chunk_size_tolerated and store either to
    the `out_file_name` if no `individual_file_prefix` is provided or
    store each individual chunk to a file starting with non-empty `individual_file_prefix`.

    Args:
        input_fasta_file: Input FASTA
        chunk_size: Size of the chunks to split into.
        chunk_size_tolerated: If more flexibility allowed, use this as the maximum size of a chunk.
        out_file_name: Output FASTA to store the chunks into if no `individual_file_prefix` is provided.
        individual_file_prefix: A file path prefix including dirs and filenames part to use as a
                first part of the chunk file name.
        agp_output_file: Output AGP file to store the map for the chunking procedure if present and non-empty.
        n_sequence_len: Length of the stretch of `N`s to split at; not slitting on `N`s if 0.
        chunk_sfx: A string to put between the original sequence name and the chunk suffix.
        append_offset_to_chunk_name: Append 0-based offset in the form of `_off_{offset}` to the chunk name.
    """

    # process input fasta
    with open_gz_file(input_fasta_file) as fasta:
        logging.info(
            f"splitting sequences from '{input_fasta_file}', chunk size {chunk_size:_}, \
                splitting on {n_sequence_len} Ns (0 -- disabled)"
        )
        # do not open a joined file if you plan to open many individual ones
        with (
            individual_file_prefix
            and nullcontext(None)
            or open(out_file_name, "wt", encoding="utf-8") as out_file_joined
        ):
            agp_lines = chunk_fasta_stream(
                fasta,
                chunk_size,
                chunk_size_tolerated,
                out_file_joined,
                individual_file_prefix,
                n_sequence_len=n_sequence_len,
                chunk_sfx=chunk_sfx,
                append_offset_to_chunk_name=append_offset_to_chunk_name,
            )

        # dump AGP
        if agp_output_file:
            with open(agp_output_file, "w", encoding="utf-8") as agp_out:
                agp_out.write("\n".join(agp_lines) + "\n")


def prepare_out_dir_for_individuals(dir_name: Path, file_part: str) -> Optional[Path]:
    """Creates `dir_name` (including upstream dirs) and returns its paths with the `file_part` appended.

    Args:
        dir_name: Directory to create.
        file_part: File part to append.

    Returns:
        `dir_name` with appended `file_part`.

    Throws:
        exception if not able to create directory.
    """
    file_prefix = None
    if dir_name:
        dir_name.mkdir(parents=True, exist_ok=True)
        file_prefix = Path(dir_name, file_part)
    return file_prefix


def main() -> None:
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
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
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
        n_sequence_len=args.n_seq,
        chunk_sfx=args.chunk_sfx,
        append_offset_to_chunk_name=args.add_offset,
    )


if __name__ == "__main__":
    main()
