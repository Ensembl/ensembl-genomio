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

"""Recombines split/chunked FASTA outputs back into a single FASTA, optionally using an AGP."""

import argparse
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
import logging
from pathlib import Path
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.io.genomio.utils.agp_utils import AgpEntry, parse_agp
from ensembl.io.genomio.utils.chunk_utils import (
    get_paths_from_manifest,
    seq_description_without_id,
    validate_regex,
)

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

_CHUNK_RE_STRING = r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$"


@dataclass(frozen=True)
class RecordLocation:
    """Location and metadata for a FASTA record discovered during header indexing."""

    path: Path
    description: str


class FastaRecordCache:
    """
    Cache for FASTA records keyed by file path.

    Only one input FASTA is loaded into memory at a time. When a requested record is in a
    different file, the cache is replaced by indexing that file's records.
    """

    def __init__(self) -> None:
        self._current_path: Path | None = None
        self._records: dict[str, SeqRecord] = {}

    def get(self, record_id: str, loc: RecordLocation) -> SeqRecord:
        """
        Return a record by ID from the file indicated by ``loc``.

        Args:
            record_id: FASTA record identifier to retrieve.
            loc: Location/metadata for the file containing the record.

        Returns:
            The requested ``SeqRecord``.

        Raises:
            KeyError: If the record ID is not found in the indexed file.
            ValueError: If the file contains duplicate record IDs.
            OSError: If the file cannot be opened/read.
        """
        if self._current_path is None or self._current_path != loc.path:
            self._load_file(loc.path)

        rec = self._records.get(record_id)
        if rec is None:
            raise KeyError(f"Record '{record_id}' not found in indexed file '{loc.path}'")
        return rec

    def _load_file(self, path: Path) -> None:
        logging.debug("Loading FASTA records from %s", path)
        records: dict[str, SeqRecord] = {}
        with open_gz_file(path) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                if rec.id in records:
                    raise ValueError(f"Duplicate record id '{rec.id}' within file '{path}'")
                records[rec.id] = rec
        self._current_path = path
        self._records = records


def _parse_fasta_headers(path: Path) -> Iterator[tuple[str, str]]:
    """
    Iterates over FASTA headers, yielding (record_id, description) tuples.

    Yields:
        Tuples of (record_id, description_without_id).
    """
    with open_gz_file(path) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            yield record.id, seq_description_without_id(record)


def _build_index(
    chunk_re: re.Pattern[str],
    fasta_paths: Iterable[Path],
) -> tuple[dict[str, RecordLocation], dict[str, int], dict[str, list[tuple[int, str]]]]:
    """
    Builds lightweight indices from FASTA headers across the manifest files.

    The function scans each FASTA file and records:
    - ``locations``: mapping of every record ID to its source file path and description.
    - ``first_seen``: stable ordering of base IDs based on first appearance across inputs
      (used to keep output order deterministic in header-driven mode).
    - ``chunks``: mapping of base ID to a list of (start, chunk_record_id) pairs extracted
      from record IDs that match ``chunk_re``.

    Args:
        chunk_re: Regex that identifies chunk records and defines named groups:
            - ``base``: original (unchunked) record ID
            - ``start``: 0-based chunk start coordinate used for ordering/contiguity checks.
        fasta_paths: Paths to FASTA files to index (may be gzipped).

    Returns:
        A tuple ``(locations, first_seen, chunks)``.

    Raises:
        ValueError: If no FASTA headers are found, or if duplicate record IDs are encountered.
        OSError: If any input file cannot be opened/read.
    """
    locations: dict[str, RecordLocation] = {}
    first_seen: dict[str, int] = {}
    chunks: dict[str, list[tuple[int, str]]] = defaultdict(list)
    order = 0

    for fasta_path in fasta_paths:
        logging.info("Indexing headers in %s", fasta_path)
        for record_id, desc in _parse_fasta_headers(fasta_path):
            if record_id in locations:
                raise ValueError(f"Duplicate FASTA record id encountered during indexing: {record_id}")
            locations[record_id] = RecordLocation(path=fasta_path, description=desc)

            m = chunk_re.match(record_id)
            base = m.group("base") if m else record_id
            if base not in first_seen:
                first_seen[base] = order
                order += 1
            if m:
                start = int(m.group("start"))
                chunks[base].append((start, record_id))

    if not locations:
        raise ValueError("No FASTA headers found in the manifest files.")
    return locations, first_seen, chunks


def _agp_component_seq(
    component_record: SeqRecord,
    comp_beg: int,
    comp_end: int,
    orientation: str,
    allow_revcomp: bool,
) -> Seq:
    """
    Extracts the component subsequence described by AGP coordinates.

    Args:
        component_record: The component sequence record (AGP ``component_id``).
        comp_beg: 1-based inclusive start on the component.
        comp_end: 1-based inclusive end on the component.
        orientation: '+' or '-' orientation from AGP.
        allow_revcomp: Whether reverse-complementing is permitted for '-' orientation.

    Returns:
        The extracted component sequence (possibly reverse-complemented).

    Raises:
        ValueError: If orientation is invalid, or if orientation is '-' when ``allow_revcomp`` is False.
    """
    sub: Seq = component_record.seq[comp_beg - 1 : comp_end]

    if orientation == "+":
        return sub
    if orientation == "-":
        if not allow_revcomp:
            raise ValueError(
                f"AGP has '-' orientation for component '{component_record.id}', "
                "but --allow-revcomp is not enabled."
            )
        return sub.reverse_complement()

    raise ValueError(f"Invalid AGP orientation '{orientation}' for component '{component_record.id}'")


def _records_from_agp(
    agp_entries: dict[str, list[AgpEntry]],
    locations: dict[str, RecordLocation],
    cache: FastaRecordCache,
    allow_revcomp: bool,
) -> Iterator[SeqRecord]:
    """
    Yields recombined records using AGP-driven reconstruction.

    For each AGP object, components are ordered by (record_start, part_number) and concatenated.
    Contiguity is enforced on the AGP object coordinates.

    Args:
        agp_entries: Mapping of AGP object ID -> list of component entries.
        locations: Mapping of record ID -> source path/description from header indexing.
        cache: Record cache used to load component sequences on demand.
        allow_revcomp: Whether reverse-complementing is permitted for '-' orientation.

    Yields:
        Reconstructed ``SeqRecord`` objects, one per AGP object.

    Raises:
        KeyError: If the AGP references a component ID not present in ``locations``.
        ValueError: If AGP parts are non-contiguous or extracted component lengths do not match AGP spans.
    """
    for record_id, parts in agp_entries.items():
        sorted_parts = sorted(parts, key=lambda p: (p.part_start, p.part_number))

        sequence_parts: list[Seq] = []
        expected_next: int | None = None
        description: str | None = None

        for part in sorted_parts:
            loc = locations.get(part.part_id)
            if loc is None:
                raise KeyError(
                    f"AGP references component_id '{part.part_id}' not found in indexed FASTA headers."
                )

            rec = cache.get(part.part_id, loc)
            sequence_part = _agp_component_seq(
                rec, part.part_start, part.part_end, part.orientation, allow_revcomp
            )
            expected_len = part.part_end - part.part_start + 1
            if len(sequence_part) != expected_len:
                raise ValueError(
                    f"Length mismatch for object '{record_id}' component '{part.part_id}': "
                    f"AGP object span={expected_len}, extracted component length={len(sequence_part)}"
                )

            if expected_next is None:
                expected_next = part.record_start
            if part.record_start != expected_next:
                raise ValueError(
                    f"Non-contiguous AGP for '{record_id}': expected next start {expected_next}, "
                    f"got {part.record_start}"
                )
            expected_next = part.record_end + 1

            sequence_parts.append(sequence_part)
            if description is None:
                description = loc.description

        assembled_sequence = Seq("").join(sequence_parts)
        yield SeqRecord(assembled_sequence, id=record_id, description=description or record_id)


def _records_from_headers(
    locations: dict[str, RecordLocation],
    first_seen: dict[str, int],
    chunks: dict[str, list[tuple[int, str]]],
    cache: FastaRecordCache,
) -> Iterator[SeqRecord]:
    """
    Yields recombined records using header-driven reconstruction.

    If a base ID has chunk records, chunk sequences are concatenated in increasing start order and
    contiguity is enforced (each chunk start must equal the previous end). If a base ID has no
    chunks, the record is passed through unchanged.

    Output order is stable based on first appearance of each base ID during indexing.

    Args:
        locations: Mapping of record ID -> source path/description.
        first_seen: Stable ordering index for each base ID.
        chunks: Mapping of base ID -> list of (start, chunk_record_id).
        cache: Record cache used to load sequences on demand.

    Yields:
        Reconstructed ``SeqRecord`` objects.

    Raises:
        KeyError: If a referenced base/chunk record is missing from ``locations``.
        ValueError: If chunk coordinates are non-contiguous.
    """
    all_bases = sorted(first_seen.keys(), key=lambda b: first_seen[b])

    for base in all_bases:
        chunk_list = chunks.get(base, [])
        if not chunk_list:
            loc = locations.get(base)
            if loc is None:
                raise KeyError(f"Base record '{base}' not found in indexed headers.")
            rec = cache.get(base, loc)
            yield SeqRecord(rec.seq, id=rec.id, description=loc.description)
            continue

        chunk_list_sorted = sorted(chunk_list, key=lambda x: x[0])

        parts: list[Seq] = []
        prev_end: int | None = None
        description: str | None = None

        for start, chunk_id in chunk_list_sorted:
            loc = locations.get(chunk_id)
            if loc is None:
                raise KeyError(f"Chunk record '{chunk_id}' not found in indexed headers.")
            rec = cache.get(chunk_id, loc)

            if prev_end is None:
                parts.append(rec.seq)
                prev_end = start + len(rec.seq)
                description = loc.description
                continue

            if start != prev_end:
                raise ValueError(
                    f"Non-contiguous chunks for '{base}': expected start {prev_end}, got {start} "
                    f"(chunk id={chunk_id})"
                )
            parts.append(rec.seq)
            prev_end = start + len(rec.seq)

        assembled = Seq("").join(parts)
        yield SeqRecord(assembled, id=base, description=description or base)


def recombine_fasta(
    fasta_manifest: Path,
    out_fasta: Path,
    chunk_re: re.Pattern[str] = re.compile(_CHUNK_RE_STRING),
    agp_file: Path | None = None,
    allow_revcomp: bool = False,
) -> None:
    """
    Recombines split/chunked FASTA files listed in a manifest into a single FASTA.

    Inputs may be plain FASTA or gzipped FASTA. Reconstruction uses one of two modes:

    1) AGP-driven (recommended):
       Assembles output sequences per AGP object using component coordinates and ordering.

    2) Header-driven (no AGP):
       Detects chunk records by matching ``chunk_re`` against record IDs, groups by ``base``, sorts by
       ``start``, enforces contiguity, and concatenates chunks. Unchunked records pass through unchanged.

    Args:
        fasta_manifest: Path to a text file containing one FASTA path per line.
        out_fasta: Destination path for the recombined FASTA output.
        chunk_re: Regex used for header-driven mode; must define named groups ``base`` and ``start``.
        agp_file: Optional AGP file. If provided, AGP-driven reconstruction is used.
        allow_revcomp: Allow reverse-complementing components when AGP orientation is '-'.

    Raises:
        ValueError: If inputs contain no records, contain duplicate IDs, or chunk/AGP coordinates are invalid.
        KeyError: If AGP references a component ID that is not present in the FASTA inputs.
        OSError: If input/output files cannot be read/written.
    """
    fasta_paths = get_paths_from_manifest(fasta_manifest)

    # Build lightweight index (headers only)
    locations, first_seen, chunks = _build_index(chunk_re, fasta_paths)

    cache = FastaRecordCache()

    if agp_file is not None:
        agp_records = parse_agp(agp_file, allow_revcomp)
        record_iter = _records_from_agp(agp_records, locations, cache, allow_revcomp)
    else:
        record_iter = _records_from_headers(locations, first_seen, chunks, cache)

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as out_fh:
        n = 0
        for rec in record_iter:
            SeqIO.write(rec, out_fh, "fasta")
            n += 1
    logging.info("Wrote %d records to %s", n, out_fasta)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "--fasta-manifest",
        metavar="TXT",
        required=True,
        help="Manifest file containing paths of FASTA files to be combined (files may be gzipped).",
    )
    parser.add_argument_dst_path(
        "--out-fasta",
        metavar="FASTA",
        required=True,
        help="Path for recombined FASTA file.",
    )
    parser.add_argument_src_path(
        "--agp-file",
        metavar="AGP",
        default=argparse.SUPPRESS,
        help="Optional AGP file; if provided, reconstruction uses AGP ordering.",
    )
    parser.add_argument(
        "--allow-revcomp",
        action="store_true",
        help="Allow reverse-complementing components if AGP orientation is '-'",
    )
    parser.add_argument(
        "--chunk-id-regex",
        type=str,
        metavar="REGEX",
        default=_CHUNK_RE_STRING,
        help=(
            "Regex used to identify chunked records and extract coordinates. "
            "Must define named groups 'base' and 'start', "
            "e.g. --chunk-id-regex '^(?P<base>.+)_(?P<start>\\d+)$'."
        ),
    )
    parser.add_log_arguments()

    args = parser.parse_args(argv)
    init_logging_with_args(args)

    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    chunk_re = validate_regex(args.chunk_id_regex)
    try:
        recombine_fasta(
            fasta_manifest=args.fasta_manifest,
            out_fasta=args.out_fasta,
            agp_file=getattr(args, "agp_file", None),
            allow_revcomp=args.allow_revcomp,
            chunk_re=chunk_re,
        )
    except Exception:
        logging.exception(f"Error recombining FASTA from files in {args.fasta_manifest}")
        raise
