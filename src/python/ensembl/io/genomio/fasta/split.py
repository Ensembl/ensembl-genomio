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

"""Split a FASTA file into multiple FASTA files, optionally chunking long sequences."""

import argparse
import logging
from pathlib import Path
import re
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


_NUMERIC_DIR_RE = re.compile(r"^[1-9]\d*$")


def _get_fasta_basename(fasta: Path) -> str:
    """Returns base name of file stripped of suffixes"""
    filename = fasta.name
    filename = filename.removesuffix(".gz")
    basename = filename.rsplit(".", 1)[0] if "." in filename else filename
    return basename


class OutputWriter:
    """
    Write split FASTA outputs and (optionally) an AGP file.

    The writer manages:
    - output directory creation/cleanup (lazy, per-directory),
    - output file naming (optionally unique across directories),
    - record and length counters used to decide when to roll over to a new file,
    - an optional AGP v2.0 file describing the mapping from original sequences
      to output contigs/chunks.

    Notes:
        Output layout is controlled by:
        - ``max_files_per_directory``: how many FASTA files to write per directory before incrementing
          the directory index.
        - ``max_dirs_per_directory``: how directory indices are expanded into a multi-level path
          (base-N style).
        - ``unique_file_names``: whether to include directory index in filenames.
    """

    def __init__(
        self,
        fasta_file: Path,
        out_dir: Path,
        write_agp: bool,
        unique_file_names: bool,
        max_files_per_directory: int | None = None,
        max_dirs_per_directory: int | None = None,
    ):
        self.basename = _get_fasta_basename(fasta_file)
        self.out_dir = out_dir
        self.agp_file = out_dir.joinpath(f"{self.basename}.agp") if write_agp else None
        self.unique_file_names = unique_file_names
        self.max_files_per_directory = max_files_per_directory
        self.max_dirs_per_directory = max_dirs_per_directory
        self.file_count = 0
        self.record_count = 0
        self.file_len = 0
        self._fh = None
        self._agp_fh = None
        self._cleaned_dirs: set[Path] = set()

        self.open_new_file()

    def _get_subdir_path(self, dir_index: int) -> Path:
        """Return the output subdirectory path for a given directory index.

        Args:
            dir_index: Zero-based directory index computed from file count.

        Returns:
            A Path under ``out_dir`` into which output files are written.
        """
        parts = []
        max_dirs = self.max_dirs_per_directory
        if max_dirs is None:
            parts.append("1")
        else:
            n = dir_index
            while n > 0:
                parts.append(str((n - 1) % max_dirs + 1))
                n = (n - 1) // max_dirs

        parts.reverse()
        return self.out_dir.joinpath(*parts)

    def _get_file_and_dir_index(self) -> tuple[int, int]:
        """Compute the file index within a directory and the directory index.

        ``file_count`` increments monotonically for each output file. If
        ``max_files_per_directory`` is set, files are grouped into directories such
        that each directory contains at most that many files.

        Returns:
            (file_index, dir_index) where:
            - file_index is 1-based within the directory, and
            - dir_index is 1-based across directories.
        """
        max_files = self.max_files_per_directory
        if max_files is None:
            return self.file_count, 0
        adjusted_count = self.file_count - 1
        return (adjusted_count % max_files + 1, adjusted_count // max_files + 1)

    def _get_path_for_next_file(self) -> Path:
        """Computes path for the next output file."""
        self.file_count += 1
        file_index, dir_index = self._get_file_and_dir_index()
        subdir_path = self._get_subdir_path(dir_index)
        subdir_path.mkdir(parents=True, exist_ok=True)

        if self.unique_file_names:
            file_name = f"{self.basename}.{dir_index}.{file_index}.fa"
        else:
            file_name = f"{self.basename}.{file_index}.fa"
        return subdir_path.joinpath(file_name)

    def add_agp_entry(
        self,
        object_id: str,
        start: int,
        end: int,
        part_nr: int,
        part_id: str,
        part_length: int,
    ) -> None:
        """
        Write a single AGP v2.0 component line for a chunk/contig.
        Coordinates written to AGP are 1-based and inclusive.

        Args:
            object_id: The original input sequence ID (AGP 'object').
            start: Start coordinate on the object (1-based, inclusive).
            end: End coordinate on the object (1-based, inclusive).
            part_nr: Component part number for this object (starts at 1 per object).
            part_id: Output contig/chunk identifier (AGP 'component_id').
            part_length: Length of the component in bases.
        """
        if self._agp_fh is None:
            return
        line = f"{object_id}\t{start}\t{end}\t{part_nr}\tW\t{part_id}\t1\t{part_length}\t+\n"
        self._agp_fh.write(line)

    def create_agp_file(self) -> None:
        """Creates the AGP file for recording sequence chunking."""
        if self.agp_file is None:
            return
        try:
            self._agp_fh = open(self.agp_file, "w")
        except OSError as e:
            raise RuntimeError(f"Failed to open AGP file '{self.agp_file}'") from e
        self._agp_fh.write("# AGP-version 2.0\n")
        logging.info(f"Created AGP file '{self.agp_file}'")

    def open_new_file(self) -> None:
        """Closes current file (if any) and opens a new output file."""
        if self._fh is not None:
            self._fh.close()

        path = self._get_path_for_next_file()
        try:
            self._fh = open(path, "w")
            logging.debug(f"Opened output file '{path}' for writing")
        except OSError as e:
            raise RuntimeError(f"Failed to open output file '{path}'") from e
        self.record_count = 0
        self.file_len = 0

    def write_record(self, record: SeqRecord) -> None:
        """Writes a SeqRecord to the current output file and update counters."""
        SeqIO.write(record, self._fh, "fasta")
        self.record_count += 1
        self.file_len += len(record.seq)

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        if self._agp_fh is not None:
            self._agp_fh.close()
            self._agp_fh = None


def _check_contents_deletable(dir: Path, output_file_re: re.Pattern[str]) -> None:
    for p in dir.rglob("*"):
        if not p.is_dir():
            if not output_file_re.match(p.name):
                msg = (
                    "Unexpected file identified amongst existing output, cleanup of existing files "
                    f"failed: {dir.parent.name}/{p.relative_to(dir.parent)}"
                )
                raise RuntimeError(msg)


def clean_previous_output(fasta_file: Path, out_dir: Path) -> None:
    """Checks for existing output and removes if no unexpected files encountered."""
    logging.info(f"Cleaning outputs from previous runs under '{out_dir}'")
    if not out_dir.exists():
        return

    # Identify top-level numeric children of out_dir (e.g. "1", "2", ...)
    top_level_dirs = [p for p in out_dir.iterdir() if p.is_dir() and _NUMERIC_DIR_RE.match(p.name)]

    basename = _get_fasta_basename(fasta_file)
    b = re.escape(basename)
    output_file_re = re.compile(rf"^{b}\.[1-9]\d*(?:\.[1-9]\d*)?\.fa$")
    for d in top_level_dirs:
        _check_contents_deletable(d, output_file_re)
    for d in top_level_dirs:
        shutil.rmtree(d)
        logging.info(f"Deleted existing output directory '{d}'.")

    agp_path = out_dir / f"{basename}.agp"
    if agp_path.exists():
        agp_path.unlink()
        logging.info(f"Deleted existing AGP file '{agp_path}'.")


def _description_without_id(record: SeqRecord) -> str:
    """Removes ID from FASTA record description"""
    desc = record.description or ""
    if desc == record.id:
        return ""
    if desc.startswith(record.id):
        return desc[len(record.id) + 1 :]
    return desc


def split_fasta(
    fasta_file: Path,
    out_dir: Path,
    write_agp: bool,
    delete_existing_files: bool,
    unique_file_names: bool,
    force_max_seq_length: bool,
    max_seqs_per_file: int | None,
    max_seq_length_per_file: int | None,
    min_chunk_length: int | None,
    max_files_per_directory: int | None = None,
    max_dirs_per_directory: int | None = None,
) -> None:
    """
    Reads an input FASTA (optionally gzipped) and writes one or more FASTA files to an output
    directory. The number of output files is determined by:

    - maximum number of records per file (``max_seqs_per_file``), and/or
    - maximum cumulative sequence length per file (``max_seq_length_per_file``).

    If ``force_max_seq_length`` is enabled, individual sequences longer than ``max_seq_length_per_file`` are
    split into chunks. When chunking, a final remainder chunk shorter than ``min_chunk_length`` can be merged
    into the previous chunk.

    Optionally, an AGP v2.0 file can be written describing how each input sequence maps to output contigs/chunks.
    """

    # Do nothing if file size is 0
    if fasta_file.stat().st_size == 0:
        logging.info(f"Input FASTA '{fasta_file}' is empty; nothing to do")
        return

    if delete_existing_files and out_dir.exists():
        clean_previous_output(fasta_file, out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    writer = OutputWriter(
        fasta_file=fasta_file,
        out_dir=out_dir,
        write_agp=write_agp,
        unique_file_names=unique_file_names,
        max_files_per_directory=max_files_per_directory,
        max_dirs_per_directory=max_dirs_per_directory,
    )

    if write_agp:
        writer.create_agp_file()

    with open_gz_file(fasta_file) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            seq_len = len(record.seq)

            if max_seqs_per_file is not None and writer.record_count >= max_seqs_per_file:
                writer.open_new_file()

            if max_seq_length_per_file is None or writer.file_len + seq_len <= max_seq_length_per_file:
                writer.write_record(record)
                if write_agp:
                    writer.add_agp_entry(record.id, 1, seq_len, 1, record.id, seq_len)
                continue

            if force_max_seq_length and seq_len > max_seq_length_per_file:
                starts = list(range(0, seq_len, max_seq_length_per_file))
                ends = [min(s + max_seq_length_per_file, seq_len) for s in starts]

                if min_chunk_length is not None and len(starts) > 1:
                    last_chunk_len = ends[-1] - starts[-1]
                    if last_chunk_len < min_chunk_length:
                        logging.warning(
                            f"Length of last chunk of record '{record.id}' is {last_chunk_len}, lower than "
                            f"min_chunk_length: {min_chunk_length}; merging with previous chunk"
                        )
                        ends[-2] = seq_len
                        starts.pop()
                        ends.pop()

                for i, (start, end) in enumerate(zip(starts, ends), start=1):
                    chunk_seq = record.seq[start:end]
                    chunk_record = SeqRecord(
                        chunk_seq,
                        id=f"{record.id}_chunk_start_{start}",
                        description=_description_without_id(record),
                    )
                    if writer.record_count > 0:
                        writer.open_new_file()
                    writer.write_record(chunk_record)

                    if write_agp:
                        writer.add_agp_entry(
                            record.id,
                            start + 1,
                            end,
                            i,
                            chunk_record.id,
                            len(chunk_seq),
                        )
            else:
                logging.warning(
                    f"Record {record.id} length {seq_len} exceeds max_seq_length_per_file "
                    f"{max_seq_length_per_file} but chunking not enabled"
                )
                if writer.record_count > 0:
                    writer.open_new_file()
                writer.write_record(record)
                if write_agp:
                    writer.add_agp_entry(record.id, 1, seq_len, 1, record.id, seq_len)
    writer.close()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--fasta-file",
        metavar="FASTA",
        required=True,
        help="Input raw or compressed FASTA file containing sequences to split.",
    )
    parser.add_argument_dst_path(
        "--out-dir",
        metavar="DIR",
        default=None,
        help="Top-level output directory (default: input FASTA directory).",
    )
    parser.add_argument(
        "--write-agp",
        action="store_true",
        help="Write AGP file describing the splits.",
    )
    parser.add_argument(
        "--delete-existing-files",
        action="store_true",
        help="Delete existing files within computed output dirs.",
    )
    parser.add_argument(
        "--unique-file-names",
        action="store_true",
        help="Make output file names unique across dirs by including dir_index.",
    )
    parser.add_argument(
        "--force-max-seq-length",
        action="store_true",
        help="Chunk single sequences longer than max-seq-length-per-file.",
    )
    parser.add_numeric_argument(
        "--max-seqs-per-file",
        type=int,
        metavar="N",
        min_value=1,
        help="Maximum number of sequences per output file.",
    )
    parser.add_numeric_argument(
        "--max-seq-length-per-file",
        type=int,
        metavar="BP",
        min_value=1,
        help="Maximum cumulative sequence length per output file.",
    )
    parser.add_numeric_argument(
        "--min-chunk-length",
        type=int,
        min_value=1,
        metavar="BP",
        help="Minimum length of a chunk allowed as a remainder.",
    )
    parser.add_numeric_argument(
        "--max-files-per-directory",
        type=int,
        min_value=1,
        metavar="N",
        help="Maximum files per directory before moving to next computed directory level.",
    )
    parser.add_numeric_argument(
        "--max-dirs-per-directory",
        type=int,
        min_value=1,
        metavar="N",
        help="Maximum subdirectories per directory level.",
    )
    parser.add_log_arguments()

    args = parser.parse_args(argv)
    if args.min_chunk_length is not None and args.max_seq_length_per_file is None:
        raise ValueError("--min-chunk-length requires --max-seq-length-per-file")

    init_logging_with_args(args)

    return args


def main(argv: list[str] | None = None) -> None:

    args = parse_args(argv)

    out_dir = args.out_dir if args.out_dir is not None else args.fasta_file.parent

    try:
        split_fasta(
            fasta_file=args.fasta_file,
            out_dir=out_dir,
            write_agp=args.write_agp,
            delete_existing_files=args.delete_existing_files,
            unique_file_names=args.unique_file_names,
            force_max_seq_length=args.force_max_seq_length,
            max_seqs_per_file=args.max_seqs_per_file,
            max_seq_length_per_file=args.max_seq_length_per_file,
            min_chunk_length=args.min_chunk_length,
            max_files_per_directory=args.max_files_per_directory,
            max_dirs_per_directory=args.max_dirs_per_directory,
        )
    except Exception:
        logging.exception(f"Error processing FASTA file {args.fasta_file}")
        raise
