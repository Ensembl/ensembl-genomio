# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Shared helpers for constructing GenomIO repeat feature JSON documents."""

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import logging
from pathlib import Path
from typing import Callable, TypeVar

T = TypeVar("T")

__all__ = [
    "Consensus",
    "file_last_modified_time",
    "format_parse_errors",
    "has_valid_parsed_coordinates",
    "parse_token",
]


@dataclass(frozen=True)
class Consensus:
    """Repeat consensus record used for feature-to-consensus linking."""

    name: str
    repeat_class: str
    repeat_type: str
    seq: str

    def sha256_key(self) -> str:
        """Return a normalized SHA256 digest for this consensus record.

        The digest is computed from the consensus name, repeat class, repeat type,
        and normalized sequence content.

        Returns:
            SHA256 hex digest for the consensus record.

        """
        norm_name = self.name.strip()
        norm_class = self.repeat_class.strip()
        norm_type = self.repeat_type.strip()
        norm_seq = "".join(self.seq.split()).upper()
        payload = f"{norm_name}\t{norm_class}\t{norm_type}\t{norm_seq}".encode()
        return hashlib.sha256(payload).hexdigest()


def parse_token(parser: Callable[[str], T], token: str, field_name: str, raw_line: str, path: Path) -> T:
    """Parse a field from tool output into the specified class type.

    Args:
        parser: A callable that takes a string and returns a value of type T.
        token: Raw token to parse.
        field_name: Field name used in error messages.
        raw_line: Original input line.
        path: Input file path.

    Returns:
        The parsed value.

    Raises:
        ValueError: If the token cannot be cast to the specified type.

    """
    try:
        return parser(token)
    except ValueError as exc:
        raise ValueError(f"Invalid {field_name!r} in {path}: token={token!r}, line={raw_line!r}") from exc


def format_parse_errors(parser_name: str, input_path: Path, errors: list[str]) -> str:
    """Format multiple parser errors into a single exception message."""
    return f"Found {len(errors)} errors while parsing {parser_name} in {input_path}:\n" + "\n".join(
        f"- {error}" for error in errors
    )


def file_last_modified_time(file_path: Path) -> str:
    """Return the last modified time of the given file."""
    return (
        datetime.fromtimestamp(
            file_path.stat().st_mtime,
            tz=timezone.utc,
        )
        .isoformat()
        .replace("+00:00", "Z")
    )


def has_valid_parsed_coordinates(
    input_path: Path,
    *,
    seq_region_start: int,
    seq_region_end: int,
    repeat_start: int,
    repeat_end: int,
    line: str,
) -> bool:
    """Validate parsed coordinate values for a feature.

    Args:
        input_path: Input file path used.
        seq_region_start: Start coordinate on the sequence region.
        seq_region_end: End coordinate on the sequence region.
        repeat_start: Start coordinate on the repeat consensus.
        repeat_end: End coordinate on the repeat consensus.
        line: Original input line for error reporting.

    Returns:
        `True` if all coordinates are valid, `False` if invalid but error not raised.

    Raises:
        ValueError: If sequence region coordinate values are invalid (i.e. negative, zero, or end < start).

    """
    if seq_region_start < 1 or seq_region_end < 1:
        raise ValueError(
            f"Invalid seq_region coordinates in {input_path}: "
            f"start={seq_region_start}, end={seq_region_end}, line={line!r}"
        )
    if seq_region_end < seq_region_start:
        raise ValueError(
            f"seq_region_end < seq_region_start in {input_path}: "
            f"start={seq_region_start}, end={seq_region_end}, line={line!r}"
        )

    no_warnings = True
    if repeat_start < 1 or repeat_end < 1:
        logging.warning(
            f"Invalid repeat coordinates in {input_path}: "
            f"repeat_start={repeat_start}, repeat_end={repeat_end}, line={line!r}"
        )
        no_warnings = False
    if repeat_end < repeat_start:
        logging.warning(
            f"repeat_end < repeat_start in {input_path}: "
            f"repeat_start={repeat_start}, repeat_end={repeat_end}, line={line!r}"
        )
        no_warnings = False

    return no_warnings
