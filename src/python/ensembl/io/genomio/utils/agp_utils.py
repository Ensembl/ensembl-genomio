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
"""Utils to deal with AGP files."""

__all__ = ["AgpEntry", "build_component_index", "lift_range", "parse_agp"]

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ensembl.utils.archive import open_gz_file


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


def build_component_index(agp_entries: dict[str, list[AgpEntry]]) -> dict[str, list[AgpEntry]]:
    """component_id -> list[AgpEntry] mapping (a component should usually map to exactly one object range)."""
    by_component: dict[str, list[AgpEntry]] = defaultdict(list)
    for _, parts in agp_entries.items():
        for p in parts:
            by_component[p.part_id].append(p)

    for comp_id, lst in by_component.items():
        by_component[comp_id] = sorted(lst, key=lambda e: (e.record, e.record_start, e.part_number))
    return by_component


def lift_range(part: AgpEntry, start: int, end: int, allow_revcomp: bool) -> tuple[str, int, int]:
    """Lift 1-based inclusive component coords -> (object_id, start, end) in object coords."""
    if start > end:
        raise ValueError(f"Range start > end: {start} > {end}")

    if start < part.part_start or end > part.part_end:
        raise ValueError(
            f"Range {start}-{end} is outside component span {part.part_start}-{part.part_end} for '{part.part_id}'"
        )

    if part.orientation == "+":
        obj_start = part.record_start + (start - part.part_start)
        obj_end = part.record_start + (end - part.part_start)
        return part.record, obj_start, obj_end

    if part.orientation == "-":
        if not allow_revcomp:
            raise ValueError(
                f"AGP has '-' orientation for component '{part.part_id}', but --allow-revcomp is not enabled."
            )
        obj_start = part.record_start + (part.part_end - end)
        obj_end = part.record_start + (part.part_end - start)
        return part.record, obj_start, obj_end

    raise ValueError(f"Invalid AGP orientation '{part.orientation}' for component '{part.part_id}'")


def parse_agp(agp_file: Path, allow_revcomp: bool) -> dict[str, list[AgpEntry]]:
    """
    Parses an AGP v2.x file into per-object component entries.

    Supported subset:
      - component type 'W' only (sequence components)
      - orientation '+' always, '-' only if `allow_revcomp=True`
    """
    agp_records: dict[str, list[AgpEntry]] = defaultdict(list)

    with open_gz_file(agp_file) as fh:
        for line in fh:
            if isinstance(line, bytes):
                line = line.decode("utf-8")
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
