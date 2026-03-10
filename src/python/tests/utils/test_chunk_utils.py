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

from contextlib import nullcontext as does_not_raise
from pathlib import Path
import pytest
from pytest import param
from typing import ContextManager


from ensembl.io.genomio.utils import chunk_utils


@pytest.mark.parametrize(
    "test_dir_name, manifest_name, expectation",
    [
        param(
            "preserves_order",
            "manifest.txt",
            does_not_raise(["b.txt", "a.txt"]),
            id="preserves_order",
        ),
        param(
            "ignores_comments_and_blank_lines",
            "manifest.txt",
            does_not_raise(["a.txt", "b.txt.gz"]),
            id="ignores_comments_and_blank_lines",
        ),
        param(
            "missing_entry",
            "manifest.txt",
            pytest.raises(FileNotFoundError, match=r"line 2"),
            id="missing_entry_raises_with_line_number",
        ),
        param(
            "directory_entry",
            "manifest.txt",
            pytest.raises(ValueError, match=r"Manifest entry.*line 1"),
            id="directory_entry_raises",
        ),
        param(
            "manifest_is_directory",
            "manifest_dir",
            pytest.raises(ValueError, match="Manifest is not a file"),
            id="manifest_must_be_a_file",
        ),
    ],
)
def test_get_paths_from_manifest(
    data_dir: Path,
    test_dir_name: str,
    manifest_name: str,
    expectation: ContextManager,
) -> None:
    """
    Tests the chunk_utils.get_paths_from_manifest() function.

    Args:
        data_dir: Module's test data directory fixture.
        test_dir_name: Name of the subdirectory within the test data directory that contains the manifest
            and its referenced files.
        manifest_name: Name of the manifest file within the test subdirectory.
        expectation: Context manager for the expected outcome of the test (exception or not).
    """
    test_dir = data_dir / test_dir_name
    manifest = test_dir / manifest_name

    with expectation as expected:
        paths = chunk_utils.get_paths_from_manifest(manifest)
        assert paths == [(test_dir / rel_path).resolve() for rel_path in expected]


def test_get_paths_from_manifest_absolute_path(tmp_path: Path) -> None:
    """
    Tests that `chunk_utils.get_paths_from_manifest()` correctly resolves absolute paths in the manifest file.

    Args:
        tmp_path: Test's unique temporary directory fixture.
    """
    a = tmp_path / "a.txt"
    a.touch()

    manifest = tmp_path / "manifest.txt"
    manifest.write_text(f"{a.resolve()}\n", encoding="utf-8")

    paths = chunk_utils.get_paths_from_manifest(manifest)
    assert paths == [a.resolve()]


@pytest.mark.parametrize(
    "pattern, text, expected_groups, expectation",
    [
        param(
            r"^(?P<base>.+)\.part\.(?P<start>\d+)_(?P<end>\d+)$",
            "contigA.part.11_20",
            {"base": "contigA", "start": "11", "end": "20"},
            does_not_raise(),
            id="accepts_alternative_pattern",
        ),
        param(
            r"^(?P<base>.+)_(?P<start>\d+$",
            None,
            None,
            pytest.raises(ValueError, match="Invalid --chunk-id-regex"),
            id="rejects_invalid_pattern",
        ),
        param(
            r"^(.+)_(\d+)$",
            None,
            None,
            pytest.raises(ValueError, match="must define named capture groups"),
            id="requires_named_groups",
        ),
        param(
            r"^(?P<base>.+)_(?P<start>\d+)$",
            "X_12",
            {"base": "X", "start": "12"},
            does_not_raise(),
            id="accepts_named_groups_pattern",
        ),
    ],
)
def test_validate_regex(
    pattern: str,
    text: str | None,
    expected_groups: dict[str, str] | None,
    expectation: ContextManager,
) -> None:
    """
    Tests the `chunk_utils.validate_regex()` function.

    Args:
        pattern: The regex pattern to validate.
        text: An optional text to match against the regex to verify the capture groups.
        expected_groups: The expected values of the 'base' and 'start' capture groups if `text` is provided.
        expectation: Context manager for the expected outcome of the test (exception or not).
    """
    with expectation:
        chunk_re = chunk_utils.validate_regex(pattern)

        if text is not None:
            match = chunk_re.match(text)
            assert match is not None
            assert expected_groups is not None
            for group, value in expected_groups.items():
                assert match.group(group) == value
