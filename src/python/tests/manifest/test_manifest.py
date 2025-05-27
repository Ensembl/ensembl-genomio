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
"""Unit testing of `ensembl.io.genomio.manifest.manifest` module."""

from contextlib import nullcontext as no_raise
import json
from pathlib import Path
from shutil import copytree
from typing import Callable, ContextManager

import pytest

from ensembl.io.genomio.manifest.manifest import Manifest, ManifestError


_CONTENT_MD5 = "45685e95985e20822fb2538a522a5ccf"


@pytest.mark.dependency
def test_init(tmp_path: Path) -> None:
    """Test `Manifest.__init__()`."""
    _ = Manifest(tmp_path)


@pytest.mark.parametrize(
    ("file_name", "expected_name"),
    [
        pytest.param("gene_models.gff3", "gff3", id="Recognised name and suffix"),
        pytest.param("foobar.gff3", "gff3", id="Recognised suffix"),
        pytest.param("genome.json", "genome", id="Recognised name"),
        pytest.param("foobar_seq_region.json", "seq_region", id="Recognised part of the name"),
        pytest.param("foobar.json", "", id="Not recognized json"),
    ],
)
@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum(tmp_path: Path, file_name: str, expected_name: str) -> None:
    """Test `Manifest.get_files_checksum()`.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        file_name: File to create for the test.
        expected_name: Manifest key expected.

    """
    expected_content = {}
    if file_name:
        with Path(tmp_path / file_name).open("w") as fh:
            fh.write("CONTENT")
        if expected_name:
            expected_content = {expected_name: {"file": file_name, "md5sum": _CONTENT_MD5}}
    manifest = Manifest(tmp_path)
    assert manifest.get_files_checksums() == expected_content
    assert manifest.files == expected_content


@pytest.mark.parametrize(
    ("file_names", "expected_content", "expected"),
    [
        pytest.param(
            ["link1.agp", "link2.agp", "link3.agp"],
            {
                "agp": {
                    "link1": {"file": "link1.agp", "md5sum": _CONTENT_MD5},
                    "link2": {"file": "link2.agp", "md5sum": _CONTENT_MD5},
                    "link3": {"file": "link3.agp", "md5sum": _CONTENT_MD5},
                },
            },
            no_raise(),
            id="3 agp files, different names",
        ),
        pytest.param(
            ["a_agp.agp", "b_agp.agp"],
            {
                "agp": {
                    "agp": {"file": "a_agp.agp", "md5sum": _CONTENT_MD5},
                    "agp.1": {"file": "b_agp.agp", "md5sum": _CONTENT_MD5},
                },
            },
            no_raise(),
            id="2 agp files with same name",
        ),
        pytest.param(
            [f"{letter}_agp.agp" for letter in "abcdefghijkl"],
            {},
            pytest.raises(ValueError, match="Too many files with same name"),
            id="Too many files with same name",
        ),
    ],
)
@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum_multifiles(
    tmp_path: Path,
    file_names: list[str],
    expected_content: dict,
    expected: ContextManager,
) -> None:
    """Test `Manifest.get_files_checksum()` with several files for the same name.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        file_names: List of files to create.
        expected_content: Expected checksum dict.
        expected: Expected exception.

    """
    for file_name in file_names:
        with Path(tmp_path / file_name).open("w") as fh:
            fh.write("CONTENT")
    manifest = Manifest(tmp_path)
    with expected:
        assert manifest.get_files_checksums() == expected_content


@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum_warning_subdir(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Test `Manifest.get_files_checksum()` with a subdir.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        caplog: Log capture fixture to check for warnings.

    """
    Path(tmp_path / "sub_dir").mkdir()
    manifest = Manifest(tmp_path)
    assert not manifest.get_files_checksums()
    assert caplog.text.strip().endswith("Can't create manifest for subdirectory")


@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum_warning_empty(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Test `Manifest.get_files_checksum()` with an empty file, deleted.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        caplog: Log capture fixture to check for warnings.

    """
    empty_file = Path(tmp_path / "empty_file.txt")
    empty_file.touch()
    assert empty_file.exists()
    manifest = Manifest(tmp_path)
    assert not manifest.get_files_checksums()
    assert "Skip and delete empty file" in caplog.text
    assert not empty_file.exists()


@pytest.mark.parametrize(
    ("files", "expected_content"),
    [
        pytest.param([], {}, id="No files"),
        pytest.param(
            ["genes.gff3"], {"gff3": {"file": "genes.gff3", "md5sum": _CONTENT_MD5}}, id="1 gff3 file"
        ),
    ],
)
@pytest.mark.dependency(depends=["test_init"])
def test_create_manifest(
    tmp_path: Path,
    assert_files: Callable,
    files: list[str],
    expected_content: dict,
) -> None:
    """Test `Manifest.create()`.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        assert_files: Function to compare if two files are equal or not.
        files: List of dummy files to create for the manifest to use.
        expected_content: Expected content of the manifest files after creation.

    """
    for file_name in files:
        with Path(tmp_path / file_name).open("w") as fh:
            fh.write("CONTENT")

    manifest = Manifest(tmp_path)
    manifest.create()
    assert manifest.file_path.exists()

    expected_manifest_file = Path(tmp_path / "expected.json")
    with expected_manifest_file.open("w") as expected_fh:
        expected_fh.write(json.dumps(expected_content, sort_keys=True, indent=4))
    assert_files(manifest.file_path, expected_manifest_file)


@pytest.mark.parametrize(
    ("files_dir", "expected_files", "expected"),
    [
        pytest.param(
            "full_data",
            {"functional_annotation", "seq_region"},
            no_raise(),
            id="OK manifest with OK files",
        ),
        pytest.param("duplicates", {"agp"}, no_raise(), id="Several files for key"),
        pytest.param("", {}, pytest.raises(ManifestError), id="No manifest to load"),
        pytest.param("missing_files", {}, pytest.raises(FileNotFoundError), id="Missing files"),
        pytest.param("invalid_checksum", {}, pytest.raises(ManifestError), id="Invalid checksum"),
    ],
)
@pytest.mark.dependency(depends=["test_init"])
def test_load(
    tmp_path: Path,
    data_dir: Path,
    files_dir: str,
    expected_files: set,
    expected: ContextManager,
) -> None:
    """Test `Manifest.load()`.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        data_dir: Module's test data directory fixture.
        files_dir: Directory where test data files are copied from.
        expected_files: Set of main files expected to be loaded.
        expected: Catch an expected exception.

    """
    # Copy the files to the tmp folder
    if files_dir:
        copytree(data_dir / files_dir, tmp_path, dirs_exist_ok=True)

    with expected:
        manifest = Manifest(tmp_path)
        manifest.load()
        assert manifest.files.keys() == expected_files
