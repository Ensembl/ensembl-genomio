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

__all__ = ["OutputWriter", "split_fasta"]

import argparse
import logging
from pathlib import Path
import re
import shutil
from typing import cast, TYPE_CHECKING

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import ensembl.io.genomio
from ensembl.io.genomio.utils.chunk_utils import seq_description_without_id

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

if TYPE_CHECKING:
    from Bio.Seq import Seq
    from io import TextIOWrapper

_NUMERIC_DIR_RE = re.compile(r"^[1-9]\d*$")


def _get_fasta_basename(fasta: Path) -> str:
    """Return the file base name stripped of suffixes."""
    filename = fasta.name
    filename = filename.removesuffix(".gz")
    return filename.rsplit(".", 1)[0]


class OutputWriter:  # pylint: disable=too-many-instance-attributes
    """Write split FASTA outputs and (optionally) an AGP file.

    The writer manages:
    - output directory creation/cleanup (lazy, per-directory),
    - output file naming (optionally unique across directories),
    - record and length counters used to decide when to roll over to a new file,
    - an optional AGP v2.0 file describing the mapping from original sequences to output contigs/chunks.

    Notes:
        Output layout is controlled by:
        - ``max_files_per_directory``: how many FASTA files to write per directory before
            incrementing the directory index.
        - ``max_dirs_per_directory``: how directory indices are expanded into a multi-level path
            (base-N style).
        - ``unique_file_names``: whether to include directory index in filenames.

    """

    def __init__(
        self,
        fasta_file: Path,
        out_dir: Path,
        *,
        write_agp: bool,
        unique_file_names: bool,
        max_files_per_directory: int | None = None,
        max_dirs_per_directory: int | None = None,
    ) -> None:
        self.basename = _get_fasta_basename(fasta_file)
        self.out_dir = out_dir
        self.write_agp = write_agp
        self.agp_file: Path | None = out_dir / f"{self.basename}.agp" if self.write_agp else None
        self.unique_file_names = unique_file_names
        self.max_files_per_directory = max_files_per_directory
        self.max_dirs_per_directory = max_dirs_per_directory
        self.file_count = 0
        self.record_count = 0
        self.file_len = 0
        self._agp_fh: TextIOWrapper | None = None
        self._cleaned_dirs: set[Path] = set()

        self._create_output_file()
        if self.write_agp:
            self._create_agp_file()

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
            return self.out_dir / "1"
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
            return (self.file_count, 0)
        adjusted_count = self.file_count - 1
        return (adjusted_count % max_files + 1, adjusted_count // max_files + 1)

    def _get_path_for_next_file(self) -> Path:
        """Compute the path for the next output file."""
        self.file_count += 1
        file_index, dir_index = self._get_file_and_dir_index()
        subdir_path = self._get_subdir_path(dir_index)
        subdir_path.mkdir(parents=True, exist_ok=True)

        if self.unique_file_names:
            file_name = f"{self.basename}.{dir_index}.{file_index}.fa"
        else:
            file_name = f"{self.basename}.{file_index}.fa"
        return subdir_path / file_name

    def _create_output_file(self) -> None:
        """Create a new output FASTA file and update internal state accordingly."""
        path = self._get_path_for_next_file()
        try:
            self._fh = path.open("w")
        except OSError as e:
            raise RuntimeError(f"Failed to open output file '{path}'") from e
        logging.debug(f"Opened output file '{path}' for writing")
        self.record_count = 0
        self.file_len = 0

    def _create_agp_file(self) -> None:
        """Create the AGP file for recording sequence chunking.

        Raises:
            ValueError: If ``write_agp`` is True but AGP file path is not set.
            RuntimeError: If the AGP file handle is not initialized.

        """
        if self.agp_file is None:
            raise ValueError("AGP file path is not set; cannot create AGP file")
        try:
            self._agp_fh = self.agp_file.open("w")
        except OSError as e:
            raise RuntimeError(f"Failed to open AGP file '{self.agp_file}'") from e
        self._agp_fh.write("# AGP-version 2.0\n")
        logging.info(f"Created AGP file '{self.agp_file}'")

    def open_new_file(self) -> None:
        """Close the current file (if any) and open a new output file."""
        self._fh.close()
        self._create_output_file()

    def write_record(
        self,
        record: SeqRecord,
        *,
        agp_object_id: str | None = None,
        agp_start: int | None = None,
        agp_end: int | None = None,
        agp_part_nr: int | None = None,
    ) -> None:
        """Write a SeqRecord to the current output file and update counters.

        If AGP writing is enabled, also writes a corresponding AGP component line describing how the
        written record maps back to the original input sequence.

        Args:
            record: Sequence record to write to the current output FASTA file.
            agp_object_id: Original (unchunked) input sequence ID to write into the AGP ``object``
                column. Required if ``write_agp`` is True.
            agp_start: Start coordinate on the AGP object (1-based, inclusive). Required if
                ``write_agp`` is True.
            agp_end: End coordinate on the AGP object (1-based, inclusive). Required if
                ``write_agp`` is True.
            agp_part_nr: Component part number for this object (starts at 1 per object). Required if
                ``write_agp`` is True.

        Raises:
            AssertionError: If ``write_agp`` is True but any of the AGP arguments are missing.
            RuntimeError: If the AGP file handle is not initialized when attempting to write an AGP line.

        """
        SeqIO.write(record, self._fh, "fasta")
        sequence = cast("Seq", record.seq)
        self.record_count += 1
        self.file_len += len(sequence)

        if self.write_agp:
            if agp_object_id is None:
                raise AssertionError("AGP object ID must be provided if writing AGP entries")
            if agp_start is None:
                raise AssertionError("AGP start must be provided if writing AGP entries")
            if agp_end is None:
                raise AssertionError("AGP end must be provided if writing AGP entries")
            if agp_part_nr is None:
                raise AssertionError("AGP part no. must be provided if writing AGP entries")
            line = (
                f"{agp_object_id}\t{agp_start}\t{agp_end}\t{agp_part_nr}\tW\t"
                f"{record.id}\t1\t{len(sequence)}\t+\n"
            )
            if self._agp_fh is None:
                raise RuntimeError("AGP file handle is not initialized")
            self._agp_fh.write(line)

    def close(self) -> None:
        """Close any open FASTA and AGP file handles."""
        self._fh.close()
        if self._agp_fh:
            self._agp_fh.close()
            self._agp_fh = None


def _check_contents_deletable(directory: Path, output_file_re: re.Pattern[str]) -> None:
    """Check that a directory contains only expected output files (recursively)."""
    for p in directory.rglob("*"):
        if not p.is_dir() and not output_file_re.match(p.name):
            msg = (
                "Unexpected file identified amongst existing output, cleanup of existing files "
                f"failed: {directory.parent.name}/{p.relative_to(directory.parent)}"
            )
            raise RuntimeError(msg)


def _clean_previous_output(fasta_file: Path, out_dir: Path) -> None:
    """Check for existing output and removes if no unexpected files encountered."""
    logging.info(f"Cleaning outputs from previous runs under '{out_dir}'")
    if not out_dir.exists():
        return

    # Identify top-level numeric children of out_dir (e.g. "1", "2", ...)
    top_level_dirs = [p for p in out_dir.iterdir() if p.is_dir() and _NUMERIC_DIR_RE.match(p.name)]

    basename = _get_fasta_basename(fasta_file)
    output_file_re = re.compile(rf"^{re.escape(basename)}\.[1-9]\d*(\.[1-9]\d*)?\.fa$")
    for d in top_level_dirs:
        _check_contents_deletable(d, output_file_re)
    for d in top_level_dirs:
        shutil.rmtree(d)
        logging.info(f"Deleted existing output directory '{d}'.")

    agp_path = out_dir / f"{basename}.agp"
    if agp_path.exists():
        agp_path.unlink()
        logging.info(f"Deleted existing AGP file '{agp_path}'.")


def split_fasta(  # noqa: PLR0912, PLR0913
    fasta_file: Path,
    out_dir: Path | None = None,
    write_agp: bool = False,
    *,
    delete_existing_files: bool = False,
    unique_file_names: bool = False,
    force_max_seq_length: bool = False,
    max_seqs_per_file: int | None = None,
    max_seq_length_per_file: int | None = None,
    min_chunk_length: int | None = None,
    max_files_per_directory: int | None = None,
    max_dirs_per_directory: int | None = None,
) -> None:
    """Read an input FASTA (optionally gzipped) and write one or more FASTA files to an output directory.

    The number of output files is determined by:
    - maximum number of records per file (``max_seqs_per_file``), and/or
    - maximum cumulative sequence length per file (``max_seq_length_per_file``).

    If ``force_max_seq_length`` is enabled, individual sequences longer than ``max_seq_length_per_file``
    are split into chunks.  Otherwise, the whole sequence will be written into a single file. When
    chunking, a final remainder chunk shorter than ``min_chunk_length`` can be merged into the
    previous chunk.

    Optionally, an AGP v2.0 file can be written describing how each input sequence maps to output
    contigs/chunks.

    Args:
        fasta_file: Input raw or compressed FASTA file containing sequences to split.
        out_dir: Top-level output directory into which split FASTA files (and optionally an AGP
            file) are written.
        write_agp: If True, write an AGP v2.0 file describing how each input sequence maps to output
            contigs/chunks.
        delete_existing_files: If True, remove outputs from previous runs under ``out_dir`` (only if
            no unexpected files are found).
        unique_file_names: If True, include directory index in output FASTA filenames to make them
            unique across nested output directories.
        force_max_seq_length: If True, sequences longer than ``max_seq_length_per_file`` are split
            into chunks.
        max_seqs_per_file: Maximum number of sequences (records) per output FASTA file. If None, no
            record-count limit is applied.
        max_seq_length_per_file: Maximum cumulative sequence length (in bp) per output FASTA file.
            If None, no cumulative-length limit is applied.
        min_chunk_length: Minimum allowed length (in bp) of the final remainder chunk when chunking.
            If the final chunk would be shorter than this, it is merged into the previous chunk. If
            None, no minimum is applied.
        max_files_per_directory: Maximum number of FASTA files per directory before moving to the
            next computed directory level. If None, all output files are written under a single directory.
        max_dirs_per_directory: Maximum number of subdirectories per directory level when expanding
            directory indices into a multi-level path (base-N style). If None, a single directory
            level is used.

    """
    out_dir = out_dir if out_dir is not None else fasta_file.parent

    # Do nothing if file size is 0
    if fasta_file.stat().st_size == 0:
        logging.info(f"Input FASTA '{fasta_file}' is empty; nothing to do")
        return

    if delete_existing_files and out_dir.exists():
        _clean_previous_output(fasta_file, out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    writer = OutputWriter(
        fasta_file=fasta_file,
        out_dir=out_dir,
        write_agp=write_agp,
        unique_file_names=unique_file_names,
        max_files_per_directory=max_files_per_directory,
        max_dirs_per_directory=max_dirs_per_directory,
    )

    try:
        with open_gz_file(fasta_file) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                seq_len = len(record.seq)

                if max_seqs_per_file is not None and writer.record_count >= max_seqs_per_file:
                    writer.open_new_file()

                if max_seq_length_per_file is None or writer.file_len + seq_len <= max_seq_length_per_file:
                    writer.write_record(
                        record, agp_object_id=record.id, agp_start=1, agp_end=seq_len, agp_part_nr=1
                    )
                    continue

                if force_max_seq_length and seq_len > max_seq_length_per_file:
                    starts = list(range(0, seq_len, max_seq_length_per_file))
                    ends = [min(s + max_seq_length_per_file, seq_len) for s in starts]

                    if min_chunk_length is not None:
                        last_chunk_len = ends[-1] - starts[-1]
                        if last_chunk_len < min_chunk_length:
                            logging.warning(
                                f"Length of last chunk of record '{record.id}' is {last_chunk_len}, "
                                f"lower than min_chunk_length {min_chunk_length}; merging with previous chunk"
                            )
                            ends[-2] = seq_len
                            starts.pop()
                            ends.pop()

                    for i, (start, end) in enumerate(zip(starts, ends, strict=True), start=1):
                        chunk_seq = record.seq[start:end]
                        chunk_record = SeqRecord(
                            chunk_seq,
                            id=f"{record.id}_chunk_start_{start}",
                            description=seq_description_without_id(record),
                        )
                        if writer.record_count > 0:
                            writer.open_new_file()
                        writer.write_record(
                            chunk_record,
                            agp_object_id=record.id,
                            agp_start=start + 1,
                            agp_end=end,
                            agp_part_nr=i,
                        )
                else:
                    logging.warning(
                        f"Record {record.id} length {seq_len} exceeds max_seq_length_per_file "
                        f"{max_seq_length_per_file} but chunking not enabled"
                    )
                    if writer.record_count > 0:
                        writer.open_new_file()
                    writer.write_record(
                        record, agp_object_id=record.id, agp_start=1, agp_end=seq_len, agp_part_nr=1
                    )
    finally:
        writer.close()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the FASTA splitting CLI.

    Args:
        argv: Optional argument vector (excluding the program name). If None, arguments are read from
            ``sys.argv`` by argparse.

    Returns:
        Parsed argparse namespace with validated options and logging configuration applied.

    Raises:
        ValueError: If ``--min-chunk-length`` is provided without ``--max-seq-length-per-file``.

    """
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
        default=argparse.SUPPRESS,
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
        help="Delete existing files within computed output directories.",
    )
    parser.add_argument(
        "--unique-file-names",
        action="store_true",
        help="Make output file names unique across nested output directories by including directory name.",
    )
    parser.add_argument(
        "--force-max-seq-length",
        action="store_true",
        help="Chunk single sequences longer than '--max-seq-length-per-file'.",
    )
    parser.add_numeric_argument(
        "--max-seqs-per-file",
        type=int,
        min_value=1,
        metavar="N",
        default=argparse.SUPPRESS,
        help="Maximum number of sequences per output file (default: no limit).",
    )
    parser.add_numeric_argument(
        "--max-seq-length-per-file",
        type=int,
        min_value=1,
        metavar="BP",
        default=argparse.SUPPRESS,
        help="Maximum cumulative sequence length per output file (default: no limit, i.e. no chunking).",
    )
    parser.add_numeric_argument(
        "--min-chunk-length",
        type=int,
        min_value=1,
        metavar="BP",
        default=argparse.SUPPRESS,
        help="Minimum length of a chunk allowed as a remainder (default: no minimum).",
    )
    parser.add_numeric_argument(
        "--max-files-per-directory",
        type=int,
        min_value=1,
        metavar="N",
        default=argparse.SUPPRESS,
        help=(
            "Maximum files per directory before moving to next computed directory level "
            "(default: all files in top-level directory)."
        ),
    )
    parser.add_numeric_argument(
        "--max-dirs-per-directory",
        type=int,
        min_value=1,
        metavar="N",
        default=argparse.SUPPRESS,
        help=(
            "Maximum subdirectories per directory level (default: all subdirectories in top-level directory)."
        ),
    )
    parser.add_log_arguments()
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)

    args = parser.parse_args(argv)
    if hasattr(args, "min_chunk_length") and not hasattr(args, "max_seq_length_per_file"):
        raise ValueError("--min-chunk-length requires --max-seq-length-per-file")

    init_logging_with_args(args)

    return args


def main(argv: list[str] | None = None) -> None:
    """Execute the FASTA splitting function."""
    args = parse_args(argv)
    try:
        split_fasta(
            fasta_file=args.fasta_file,
            out_dir=getattr(args, "out_dir", None),
            write_agp=args.write_agp,
            delete_existing_files=args.delete_existing_files,
            unique_file_names=args.unique_file_names,
            force_max_seq_length=args.force_max_seq_length,
            max_seqs_per_file=getattr(args, "max_seqs_per_file", None),
            max_seq_length_per_file=getattr(args, "max_seq_length_per_file", None),
            min_chunk_length=getattr(args, "min_chunk_length", None),
            max_files_per_directory=getattr(args, "max_files_per_directory", None),
            max_dirs_per_directory=getattr(args, "max_dirs_per_directory", None),
        )
    except Exception:
        logging.exception(f"Error processing FASTA file {args.fasta_file}")
        raise
