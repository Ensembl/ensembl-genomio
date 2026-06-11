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

from ensembl.utils import StrPath


def get_paths_from_manifest(manifest: StrPath) -> list[Path]:
    """Parses a manifest file and return a list of validated input paths.

    Args:
        manifest: Path to a manifest file containing one input path per line.

    Returns:
        A list of resolved ``Path`` objects for each file listed in the manifest,
        preserving the manifest order.

    Raises:
        FileNotFoundError: If the manifest file or any listed entry does not exist.
        ValueError: If the manifest path is not a file or any listed entry is not a file.

    """
    manifest = Path(manifest)
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
    """Returns the FASTA record description with the record ID removed.

    Args:
        record: A ``SeqRecord`` whose ``description`` may start with ``id``.

    Returns:
        The description text with the leading identifier and following space
        removed. Returns an empty string if the description equals the id.

    """
    desc = record.description
    if record.id is None or desc == record.id:
        return ""
    return desc.removeprefix(record.id + " ")


def validate_regex(chunk_regex: str) -> re.Pattern[str]:
    """Compiles and validates a chunk-id regex.

    The compiled regex must define named capture groups ``base`` and ``start``.

    Args:
        chunk_regex: Regular expression string to compile.

    Returns:
        A compiled ``re.Pattern`` instance.

    Raises:
        ValueError: If the regex is invalid or does not contain the required
            named capture groups ``base`` and ``start``.

    """
    try:
        chunk_re = re.compile(chunk_regex)
    except re.error as e:
        raise ValueError(f"Invalid regex: {e}") from e

    if "base" not in chunk_re.groupindex or "start" not in chunk_re.groupindex:
        raise ValueError("Chunk ID regex must define named capture groups 'base' and 'start'")

    return chunk_re
