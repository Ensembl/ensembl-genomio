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
"""Unit testing of `ensembl.io.genomio.utils.chunk_utils` module."""

import pytest

from ensembl.io.genomio.utils import chunk_utils
from .._helpers import write_manifest


def test_get_paths_from_manifest_preserves_order(tmp_path):
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.touch()
    b.touch()

    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            str(b),
            str(a),
        ],
    )

    paths = chunk_utils.get_paths_from_manifest(manifest)
    assert paths == [b.resolve(), a.resolve()]


def test_get_paths_from_manifest_ignores_comments_and_blank_lines(tmp_path):
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt.gz"
    a.touch()
    b.touch()

    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            "# comment",
            "",
            str(a),
            "   ",
            str(b),
            "# trailing comment",
        ],
    )

    paths = chunk_utils.get_paths_from_manifest(manifest)
    assert paths == [a.resolve(), b.resolve()]


def test_get_paths_from_manifest_relative_paths_are_anchored_to_manifest_dir(tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    a = tmp_path / "sub" / "a.txt"
    a.touch()

    manifest = write_manifest(
        sub / "manifest.txt",
        [
            "a.txt",  # relative to manifest dir
        ],
    )

    paths = chunk_utils.get_paths_from_manifest(manifest)
    assert paths == [a.resolve()]


def test_get_paths_from_manifest_missing_entry_raises_with_line_number(tmp_path):
    a = tmp_path / "a.txt"
    a.touch()

    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            str(a),
            "missing.txt",
        ],
    )

    with pytest.raises(FileNotFoundError, match=r"line 2"):
        chunk_utils.get_paths_from_manifest(manifest)


def test_get_paths_from_manifest_directory_entry_raises(tmp_path):
    d = tmp_path / "adir"
    d.mkdir()

    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            str(d),
        ],
    )

    with pytest.raises(ValueError, match=r"Manifest entry.*line 1"):
        chunk_utils.get_paths_from_manifest(manifest)


def test_get_paths_from_manifest_manifest_must_be_a_file(tmp_path):
    manifest_dir = tmp_path / "manifest_dir"
    manifest_dir.mkdir()

    with pytest.raises(ValueError, match="Manifest is not a file"):
        chunk_utils.get_paths_from_manifest(manifest_dir)


def test_validate_regex_accepts_alternative_pattern():
    # Alternative chunk naming scheme: "<base>.part.<start>_<end>"
    chunk_re = chunk_utils.validate_regex(r"^(?P<base>.+)\.part\.(?P<start>\d+)_(?P<end>\d+)$")
    m = chunk_re.match("contigA.part.11_20")
    assert m is not None
    assert m.group("base") == "contigA"
    assert m.group("start") == "11"


def test_validate_regex_rejects_invalid_pattern():
    with pytest.raises(ValueError, match="Invalid --chunk-id-regex"):
        chunk_utils.validate_regex(r"^(?P<base>.+)_(?P<start>\d+$")  # missing ')'


def test_validate_regex_requires_named_groups():
    with pytest.raises(ValueError, match="must define named capture groups"):
        chunk_utils.validate_regex(r"^(.+)_(\d+)$")

    ok = chunk_utils.validate_regex(r"^(?P<base>.+)_(?P<start>\d+)$")
    assert ok.match("X_12")
