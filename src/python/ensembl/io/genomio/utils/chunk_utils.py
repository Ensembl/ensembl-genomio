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
"""Utils to deal with processing of chunked files."""

__all__ = ["get_paths_from_manifest", "seq_description_without_id", "validate_regex"]

from pathlib import Path
import re

from Bio.SeqRecord import SeqRecord


def get_paths_from_manifest(manifest: Path) -> list[Path]:
    """
    Parses manifest file to return list of validated input paths
    """
    manifest = manifest.expanduser().resolve(strict=True)

    if not manifest.is_file():
        raise ValueError(f"Manifest is not a file: {manifest}")

    base_dir = manifest.parent
    paths: list[Path] = []

    with manifest.open() as fh:
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


def seq_description_without_id(record: SeqRecord) -> str:
    """Removes ID from FASTA record description"""
    desc = record.description
    if desc == record.id:
        return ""
    if desc.startswith(record.id):
        return desc[len(record.id) + 1 :]
    return desc


def validate_regex(chunk_regex) -> re.Pattern[str]:
    """Compiles and validates the chunk-id regex, ensuring 'base' and 'start' capture groups present."""
    try:
        chunk_re = re.compile(chunk_regex)
    except re.error as e:
        raise ValueError(f"Invalid --chunk-id-regex: {e}") from e

    if "base" not in chunk_re.groupindex or "start" not in chunk_re.groupindex:
        raise ValueError("--chunk-id-regex must define named capture groups 'base' and 'start'")

    return chunk_re
