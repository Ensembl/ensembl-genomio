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


MIN_AGP_COLUMNS = 9
SUPPORTED_COMPONENT_TYPES = {"W"}

@dataclass(frozen=True)
class AgpEntry:
    """Single AGP component entry used for coordinate liftover."""

    record: str
    record_start: int
    record_end: int
    part_number: int
    part_id: str
    part_start: int
    part_end: int
    orientation: str


def build_component_index(agp_entries: dict[str, list[AgpEntry]]) -> dict[str, list[AgpEntry]]:
    """Build an index of component IDs to AGP entries.

    Args:
        agp_entries: Mapping of object IDs to lists of AGP entries.

    Returns:
        A mapping from component ID to sorted lists of AGP entries.

    """
    component_index: dict[str, list[AgpEntry]] = defaultdict(list)
    for parts in agp_entries.values():
        for p in parts:
            component_index[p.part_id].append(p)

    for key, values in component_index.items():
        component_index[key] = sorted(
            values,
            key=lambda e: (e.record, e.record_start, e.part_number),
        )
    return component_index


def lift_range(part: AgpEntry, start: int, end: int, allow_revcomp: bool) -> tuple[str, int, int]:
    """Lift component coordinates into object coordinates.

    Args:
        part: AGP entry describing the component span.
        start: Start coordinate on the component (1-based inclusive).
        end: End coordinate on the component (1-based inclusive).
        allow_revcomp: Whether reverse-complement orientation is permitted.

    Returns:
        A tuple of object ID, object start, and object end.

    Raises:
        ValueError: If the requested range is invalid or not covered by the component.

    """
    if start > end:
        raise ValueError(f"Range start > end: {start} > {end}")

    if start < part.part_start or end > part.part_end:
        raise ValueError(
            f"Range {start}-{end} is outside component span "
            f"{part.part_start}-{part.part_end} for '{part.part_id}'"
        )

    if part.orientation == "+":
        obj_start = part.record_start + (start - part.part_start)
        obj_end = part.record_start + (end - part.part_start)
        return part.record, obj_start, obj_end

    if part.orientation == "-":
        if not allow_revcomp:
            raise ValueError(
                f"AGP contains '-' orientation for component '{part.part_id}' but processing of "
                "reverse complement AGP entries is not enabled."
            )
        obj_start = part.record_start + (part.part_end - end)
        obj_end = part.record_start + (part.part_end - start)
        return part.record, obj_start, obj_end

    raise ValueError(f"Invalid AGP orientation '{part.orientation}' for component '{part.part_id}'")


def parse_agp(agp_file: Path, allow_revcomp: bool) -> dict[str, list[AgpEntry]]:
    """Parse an AGP v2.x file into per-object component entries.

    Args:
        agp_file: Path to the input AGP file (optionally gzipped).
        allow_revcomp: Whether '-' orientations are permitted.

    Returns:
        Mapping from object ID to lists of AGP entries.

    Raises:
        ValueError: If the AGP file is invalid, unsupported, or contains no component lines.
                    All validation errors are collected and raised together.

    """
    agp_records: dict[str, list[AgpEntry]] = defaultdict(list)
    errors: list[str] = []

    with open_gz_file(agp_file) as fh:
        for line_nr, raw_line in enumerate(fh, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < MIN_AGP_COLUMNS:
                errors.append(
                    f"Line {line_nr}: Invalid AGP line (expected >= {MIN_AGP_COLUMNS} columns): {line}"
                )
                continue

            if cols[4] not in SUPPORTED_COMPONENT_TYPES:
                errors.append(f"Line {line_nr}: Unsupported AGP component type '{cols[4]}' in line: {line}")
                continue

            if not allow_revcomp and cols[8] != "+":
                errors.append(
                    f"Line {line_nr}: AGP contains '-' orientation for component '{cols[5]}' "
                    "but processing of reverse complement AGP entries is not enabled."
                )
                continue

            try:
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
            except ValueError as e:
                errors.append(f"Line {line_nr}: {e!s}")

    if errors:
        error_msg = f"AGP file '{agp_file}' had {len(errors)} error(s):\n  - " + "\n  - ".join(errors)
        raise ValueError(error_msg)

    if not agp_records:
        raise ValueError(f"AGP file '{agp_file}' contained no component lines.")
    return agp_records
