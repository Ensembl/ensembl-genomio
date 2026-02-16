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

"""Recombine split/chunked FASTA outputs back into a single FASTA, optionally using an AGP."""

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

from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


@dataclass(frozen=True)
class RecordLocation:
    path: Path
    description: str


@dataclass(frozen=True)
class AgpEntry:
    record: str
    record_start: int
    record_end: int
    part_number: int
    part_id: str
    part_start: int
    part_end: int
    orientation: str


class FastaRecordCache:
    """Cache for FASTA records keyed by file path, to keep single file's records in memory at once."""

    def __init__(self) -> None:
        self._current_path: Path | None = None
        self._records: dict[str, SeqRecord] = {}

    def get(self, record_id: str, loc: RecordLocation) -> SeqRecord:
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


def _get_fasta_paths(fasta_manifest: Path) -> list[Path]:
    """
    Parses manifest file to return list of validated input FASTA paths
    """
    fasta_manifest = fasta_manifest.expanduser().resolve(strict=True)

    if not fasta_manifest.is_file():
        raise ValueError(f"Manifest is not a file: {fasta_manifest}")

    base_dir = fasta_manifest.parent
    paths: list[Path] = []

    with fasta_manifest.open() as fh:
        for line_nr, line in enumerate(fh, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            p = Path(line).expanduser()
            if not p.is_absolute():
                p = base_dir / p
            try:
                p = p.resolve(strict=True)
            except FileNotFoundError as e:
                raise FileNotFoundError(f"Manifest entry does not exist (line {line_nr}): {line}") from e
            if not p.is_file():
                raise ValueError(f"Manifest entry is not a file (line {line_nr}): {p}")

            paths.append(p)

    return paths


def _description_without_id(record: SeqRecord) -> str:
    """Removes ID from FASTA record description"""
    if record.description.startswith(record.id):
        return record.description[len(record.id) :].lstrip()
    return record.description


def _parse_fasta_headers(path: Path) -> Iterator[tuple[str, str]]:
    """Iterate over FASTA headers, yielding (record_id, description) tuples."""
    with open_gz_file(path) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            yield record.id, _description_without_id(record)


def _build_index(
    chunk_re: re.Pattern[str],
    fasta_paths: Iterable[Path],
) -> tuple[dict[str, RecordLocation], dict[str, int], dict[str, list[tuple[int, str]]]]:
    """
    Builds an index from record_id -> (file_path, description) using headers only.

    Returns:
      - locations: record_id -> RecordLocation(path, description)
      - first_seen: base_id -> first-seen order index (stable output ordering in header mode)
      - chunks: base_id -> list of (start, chunk_record_id) for header-driven reconstruction
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


def _parse_agp(agp_file: Path, allow_revcomp: bool) -> dict[str, list[AgpEntry]]:
    """
    Parses an AGP v2.x file into per-object component entries.

    Supported subset:
      - component type 'W' only (sequence components)
      - orientation '+' always, '-' only if `allow_revcomp=True`
    """
    agp_records: dict[str, list[AgpEntry]] = defaultdict(list)

    with open_gz_file(agp_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                raise ValueError(f"Invalid AGP line (expected >= 9 columns): {line}")

            if cols[4] != "W":
                raise ValueError(f"Unsupported AGP component type '{cols[4]}' in line: {line}")

            if not allow_revcomp and cols[8] != "+":
                raise ValueError(
                    f"AGP contains '-' orientation for component '{cols[5]}' but --allow-revcomp is not enabled."
                )

            agp_records[cols[0]].append(
                AgpEntry(
                    record=cols[0],
                    record_start=int(cols[1]),
                    record_end=int(cols[2]),
                    part_number=int(cols[3]),
                    part_id=cols[5],
                    part_start=int(cols[6]),
                    part_end=int(cols[7]),
                    orientation=cols[8],
                )
            )

    if not agp_records:
        raise ValueError(f"AGP file '{agp_file}' contained no component lines.")
    return agp_records


def _agp_component_seq(
    component_record: SeqRecord,
    comp_beg: int,
    comp_end: int,
    orientation: str,
    allow_revcomp: bool,
) -> Seq:
    """
    Extracts the component subsequence described by AGP coordinates.

    Converts 1-based AGP coordinates to Python slicing and returns the subsequence.
    If orientation is '-', reverse-complement is applied when allowed.
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
    """AGP-driven reconstruction of sequences."""
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
                    f"Non-contiguous AGP for '{record_id}': expected next start {expected_next}, got {part.record_start}"
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
    Header-driven reconstruction of sequences:
      - If a base id has chunk records, concatenate in start order (enforcing contiguity).
      - Otherwise, pass through the unchunked record.
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
    chunk_re: re.Pattern[str],
    agp_file: Path | None = None,
    allow_revcomp: bool = False,
) -> None:
    """
    Recombines split/chunked FASTA files listed in manifest into a single FASTA.  Inputs may be gzipped.

    Reconstruction modes:

    1) AGP-driven (recommended)
    - Output sequences are assembled per AGP 'object' in coordinate order.

    2) Header-driven (no AGP)
    - Chunk records detected by parsing the header
    - Chunks grouped by <orig_id>, sorted by start, and concatenated.
    - Unchunked records pass through unchanged.
    """
    fasta_paths = _get_fasta_paths(fasta_manifest)

    # Build lightweight index (headers only)
    locations, first_seen, chunks = _build_index(chunk_re, fasta_paths)

    cache = FastaRecordCache()

    if agp_file is not None:
        agp_records = _parse_agp(agp_file, allow_revcomp)
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


def _validate_regex(chunk_regex) -> re.Pattern[str]:
    """Compiles and validates the chunk-id regex, ensuring 'base' and 'start' capture groups present."""
    try:
        chunk_re = re.compile(chunk_regex)
    except re.error as e:
        raise ValueError(f"Invalid --chunk-id-regex: {e}") from e

    if "base" not in chunk_re.groupindex or "start" not in chunk_re.groupindex:
        raise ValueError("--chunk-id-regex must define named capture groups 'base' and 'start'")

    return chunk_re


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
        default=None,
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
        default=r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$",
        help=(
            "Regex used to identify chunked records and extract coordinates. "
            "Must define named groups 'base' and 'start', e.g. --chunk-id-regex '^(?P<base>.+)_(?P<start>\\d+)$'."
        ),
    )

    args = parser.parse_args(argv)
    init_logging_with_args(args)

    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    chunk_re = _validate_regex(args.chunk_id_regex)
    try:
        recombine_fasta(
            fasta_manifest=args.fasta_manifest,
            out_fasta=args.out_fasta,
            agp_file=args.agp_file,
            allow_revcomp=args.allow_revcomp,
            chunk_re=chunk_re,
        )
    except Exception:
        logging.exception(f"Error recombining FASTA from files in {args.fasta_manifest}")
        raise
