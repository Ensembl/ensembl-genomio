#!/usr/bin/env python3
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

"""
Recombine split FASTA outputs back into a single FASTA, optionally using an AGP.

Reconstruction modes:

1) AGP-driven (recommended)
   - Provide --agp-file (e.g. <basename>.agp created by the splitter with --write-agp)
   - Output sequences are assembled per AGP 'object' in coordinate order.

2) Header-driven (no AGP)
   - Detect chunk records created by the splitter:
       <orig_id>_chunk_start_<0-based-start>
   - Group by <orig_id>, sort by start, and concatenate.
   - Unchunked records pass through unchanged.

Input FASTA records can be spread across multiple files and nested directories under --in-dir.
"""

import inspect
import logging
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Pattern, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.utils.archive import open_gz_file  # type: ignore
from ensembl.utils.argparse import ArgumentParser  # type: ignore
from ensembl.utils.logging import init_logging_with_args  # type: ignore

_num_re = re.compile(r"(\d+)")


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


class Params:
    """Validated configuration for recombining FASTA splits."""

    def __init__(
        self,
        in_dir: Path,
        out_fasta: Path,
        agp_file: Optional[Path] = None,
        extra_suffix: Optional[List[str]] = None,
        allow_revcomp: bool = False,
        chunk_id_regex: Optional[str] = r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$",
    ):
        self.in_dir = in_dir
        self.out_fasta = out_fasta
        self.agp_file = agp_file
        self.extra_suffixes = extra_suffix or []
        self.allow_revcomp = allow_revcomp
        self._validate_params(chunk_id_regex)

    def _validate_params(self, chunk_regex) -> None:
        try:
            self.chunk_re: Pattern[str] = re.compile(chunk_regex)

        except re.error as e:
            raise ValueError(f"Invalid --chunk-regex: {e}") from e
        if "base" not in self.chunk_re.groupindex or "start" not in self.chunk_re.groupindex:
            raise ValueError("--chunk-regex must define named capture groups 'base' and 'start'")
        if not self.in_dir.exists() or not self.in_dir.is_dir():
            raise ValueError(f"--in-dir does not exist or is not a directory: {self.in_dir}")
        if self.agp_file is not None and not self.agp_file.exists():
            raise ValueError(f"--agp-file does not exist: {self.agp_file}")


class FastaRecordCache:
    """
    Cache for FASTA records keyed by file path in order to keep at most
    one file's records in memory at once.
    """

    def __init__(self) -> None:
        self._current_path: Optional[Path] = None
        self._records: Dict[str, SeqRecord] = {}

    def get(self, record_id: str, loc: RecordLocation) -> SeqRecord:
        if self._current_path is None or self._current_path != loc.path:
            self._load_file(loc.path)

        rec = self._records.get(record_id)
        if rec is None:
            raise KeyError(f"Record '{record_id}' not found in indexed file '{loc.path}'")
        return rec

    def _load_file(self, path: Path) -> None:
        logging.debug("Loading FASTA records from %s", path)
        records: Dict[str, SeqRecord] = {}
        with open_gz_file(path) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                if rec.id in records:
                    raise ValueError(f"Duplicate record id '{rec.id}' within file '{path}'")
                records[rec.id] = rec
        self._current_path = path
        self._records = records


def _strip_fasta_suffix(name: str, suffixes: List[str]) -> str:
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _numeric_path_key(path: Path, suffixes: List[str]) -> List[str]:
    key = []
    parts = path.parts

    for i, part in enumerate(parts):
        text = part

        if i == len(parts) - 1:
            text = _strip_fasta_suffix(text, suffixes)

        for token in _num_re.split(text):
            if not token:
                continue
            if token.isdigit():
                key.append((0, int(token)))  # numbers first, numeric compare
            else:
                key.append((1, token))  # then text, lexicographic compare
    return key


def _get_fasta_paths(params: Params) -> List[Path]:
    suffixes = ["fa", "fasta", "fna"] + list(params.extra_suffixes)

    seen: Dict[Path, None] = {}
    for suffix in suffixes:
        compression_suffixes = ("",) if suffix.endswith(".gz") else ("", ".gz")
        for compression_suffix in compression_suffixes:
            pattern = f"**/*{suffix}{compression_suffix}"
            for result in params.in_dir.glob(pattern):
                if result.is_file():
                    seen[result.resolve()] = None
    return sorted(seen.keys(), key=lambda p: _numeric_path_key(p, suffixes))


def _parse_fasta_headers(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Iterate over FASTA headers, yielding (record_id, description) tuples.
    """
    with open_gz_file(path) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            yield record.id, record.description


def _build_index(
    params: Params,
    fasta_paths: Iterable[Path],
) -> Tuple[Dict[str, RecordLocation], Dict[str, int], Dict[str, List[Tuple[int, str]]]]:
    """
    Build an index from record_id -> (file_path, description) using headers only.

    Returns:
      - locations: record_id -> RecordLocation(path, description)
      - first_seen: base_id -> first-seen order index (stable output ordering in header mode)
      - chunks: base_id -> list of (start, chunk_record_id) for header-driven reconstruction
    """
    locations: Dict[str, RecordLocation] = {}
    first_seen: Dict[str, int] = {}
    chunks: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    order = 0

    for fasta_path in fasta_paths:
        logging.info("Indexing headers in %s", fasta_path)
        for record_id, desc in _parse_fasta_headers(fasta_path):
            if record_id in locations:
                raise ValueError(f"Duplicate FASTA record id encountered during indexing: {record_id}")
            locations[record_id] = RecordLocation(path=fasta_path, description=desc)

            m = params.chunk_re.match(record_id)
            base = m.group("base") if m else record_id
            if base not in first_seen:
                first_seen[base] = order
                order += 1
            if m:
                start = int(m.group("start"))
                chunks[base].append((start, record_id))

    if not locations:
        raise ValueError("No FASTA headers found under the input directory.")
    return locations, first_seen, chunks


def _parse_agp(params: Params) -> Dict[str, List[AgpEntry]]:
    agp_entries: Dict[str, List[AgpEntry]] = defaultdict(list)

    with open_gz_file(params.agp_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                raise ValueError(f"Invalid AGP line (expected >= 9 columns): {line}")

            if cols[4] != "W":
                raise ValueError(f"Unsupported AGP component type '{cols[4]}' in line: {line}")

            if not params.allow_revcomp and cols[8] != "+":
                raise ValueError(
                    f"AGP contains '-' orientation for component '{cols[5]}' but --allow-revcomp is not enabled."
                )

            agp_entries[cols[0]].append(
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

    if not agp_entries:
        raise ValueError(f"AGP file '{params.agp_file}' contained no component lines.")
    return agp_entries


def _agp_component_seq(
    component_record: SeqRecord,
    comp_beg: int,
    comp_end: int,
    orientation: str,
    allow_revcomp: bool,
) -> Seq:
    # AGP coords: 1-based inclusive -> Python slice [beg-1 : end]
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
    agp_entries: Dict[str, List[AgpEntry]],
    locations: Dict[str, RecordLocation],
    cache: FastaRecordCache,
    allow_revcomp: bool,
) -> Iterator[SeqRecord]:
    for record_id, parts in agp_entries.items():
        sorted_parts = sorted(parts, key=lambda p: (p.part_start, p.part_number))

        sequence_parts: List[Seq] = []
        expected_next: Optional[int] = None
        description: Optional[str] = None

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
                    f"Non-contiguous AGP for '{record_id}': expected next start {expected_next}, got {c.record_start}"
                )
            expected_next = part.record_end + 1

            sequence_parts.append(sequence_part)
            if description is None:
                description = loc.description

        assembled_sequence = Seq("").join(sequence_parts)
        yield SeqRecord(assembled_sequence, id=record_id, description=description or record_id)


def _records_from_headers(
    locations: Dict[str, RecordLocation],
    first_seen: Dict[str, int],
    chunks: Dict[str, List[Tuple[int, str]]],
    cache: FastaRecordCache,
) -> Iterator[SeqRecord]:
    """
    Header-driven reconstruction:
      - If a base id has chunk records, concatenate in start order (enforcing contiguity).
      - Otherwise, pass through the unchunked record.
    """
    all_bases = sorted(first_seen.keys(), key=lambda b: first_seen[b])

    for base in all_bases:
        chunk_list = chunks.get(base, [])
        if not chunk_list:
            # unchunked: base record id must exist
            loc = locations.get(base)
            if loc is None:
                raise KeyError(f"Base record '{base}' not found in indexed headers.")
            rec = cache.get(base, loc)
            yield SeqRecord(rec.seq, id=rec.id, description=loc.description)
            continue

        chunk_list_sorted = sorted(chunk_list, key=lambda x: x[0])

        parts: List[Seq] = []
        prev_end: Optional[int] = None
        description: Optional[str] = None

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


def recombine_fasta(params: Params) -> None:
    fasta_paths = _get_fasta_paths(params)
    if not fasta_paths:
        raise ValueError(f"No FASTA files found under: {params.in_dir}")

    # Build lightweight index (headers only)
    locations, first_seen, chunks = _build_index(params, fasta_paths)

    cache = FastaRecordCache()

    if params.agp_file is not None:
        agp_by_obj = _parse_agp(params)
        record_iter = _records_from_agp(agp_by_obj, locations, cache, params.allow_revcomp)
    else:
        record_iter = _records_from_headers(locations, first_seen, chunks, cache)

    params.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(params.out_fasta, "w") as out_fh:
        n = 0
        for rec in record_iter:
            # Write per-record to avoid materialising the whole iterator
            SeqIO.write(rec, out_fh, "fasta")
            n += 1
    logging.info("Wrote %d records to %s", n, params.out_fasta)


def _get_param_defaults() -> dict:
    """
    Return default values from the ``Params`` constructor signature.

    Keeps CLI help text in sync with the defaults defined in ``Params.__init__``.
    """
    signature = inspect.signature(Params.__init__)
    defaults = {}
    for name, param in signature.parameters.items():
        if name != "self" and param.default is not inspect.Parameter.empty:
            defaults[name] = param.default
    return defaults


def parse_args(argv: Optional[List[str]] = None) -> Params:
    """Parse CLI arguments and return a validated Params object."""
    defaults = _get_param_defaults()
    parser = ArgumentParser(description="Recombine split FASTA outputs into a single FASTA.")
    parser.add_argument(
        "--in-dir",
        type=Path,
        metavar="DIR",
        required=True,
        help="Directory containing split FASTA outputs (searched recursively, split files may be gzipped)",
    )
    parser.add_argument(
        "--out-fasta",
        type=Path,
        metavar="FASTA",
        required=True,
        help="Output recombined FASTA file",
    )
    parser.add_argument(
        "--agp-file",
        type=Path,
        metavar="AGP",
        help=f"Optional AGP file; if provided, reconstruction uses AGP ordering (default: {defaults['agp_file']})",
    )
    parser.add_argument(
        "--extra-suffix",
        action="append",
        metavar="SUFFIX",
        help=(
            "Additional file suffixes to search for under --in-dir (repeatable), e.g. 'fsa'. "
            "(fa/fasta/fna are searched for by default, .gz doesn't need to be specified)"
        ),
    )
    parser.add_argument(
        "--allow-revcomp",
        action="store_true",
        help=f"Allow reverse-complementing components if AGP orientation is '-' (default: {defaults['allow_revcomp']})",
    )
    parser.add_argument(
        "--chunk-id-regex",
        type=str,
        metavar="REGEX",
        help=(
            "Regex used to identify chunked records and extract coordinates. "
            "Must define named groups 'base' and 'start'. "
            f"(default: {defaults['chunk_id_regex']})"
        ),
    )

    args = parser.parse_args(argv)
    init_logging_with_args(args)

    return Params(
        in_dir=args.in_dir,
        out_fasta=args.out_fasta,
        agp_file=args.agp_file,
        extra_suffix=args.extra_suffix,
        allow_revcomp=args.allow_revcomp,
        chunk_id_regex=args.chunk_id_regex,
    )


def main(argv: Optional[List[str]] = None) -> None:
    try:
        params = parse_args(argv)
        recombine_fasta(params)
    except Exception:
        logging.exception("Error recombining FASTA")
        raise


if __name__ == "__main__":
    main()
