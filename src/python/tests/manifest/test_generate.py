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
"""Unit testing of `ensembl.io.genomio.manifest.generate` module."""

import json
from pathlib import Path
from typing import Callable

import pytest
from pytest import LogCaptureFixture, param

from ensembl.io.genomio.manifest.generate import ManifestMaker


_CONTENT_MD5 = "45685e95985e20822fb2538a522a5ccf"


@pytest.mark.dependency()
def test_init(tmp_path: Path) -> None:
    """Tests the `ManifestMaker.__init__()` method."""
    _ = ManifestMaker(tmp_path)


@pytest.mark.parametrize(
    "file_name, expected_name",
    [
        param("gene_models.gff3", "gff3", id="Recognised name and suffix"),
        param("foobar.gff3", "gff3", id="Recognised suffix"),
        param("genome.json", "genome", id="Recognised name"),
        param("foobar_seq_region.json", "seq_region", id="Recognised part of the name"),
        param("foobar.json", "", id="Not recognized json"),
    ],
)
@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum(tmp_path: Path, file_name: str, expected_name: str) -> None:
    """Tests the `ManifestMaker.get_files_checksum()` method.

    Args:
        tmp_path: Test tmp dir.
        expected_name: Manifest key expected.
        file_name: File to create for the test.

    """
    expected_content = {}
    if file_name:
        with Path(tmp_path / file_name).open("w") as fh:
            fh.write("CONTENT")
        if expected_name:
            expected_content = {expected_name: {"file": file_name, "md5sum": _CONTENT_MD5}}
    maker = ManifestMaker(tmp_path)
    test_files = maker.get_files_checksums()
    assert test_files == expected_content


@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum_multifiles(tmp_path: Path) -> None:
    """Tests the `ManifestMaker.get_files_checksum()` method with several files for the same name."""
    expected_content = {}
    files = ["link1.agp", "link2.agp", "link3.agp"]
    for file_name in files:
        with Path(tmp_path / file_name).open("w") as fh:
            fh.write("CONTENT")
    expected_content = {
        "agp": {
            "link1": {"file": "link1.agp", "md5sum": _CONTENT_MD5},
            "link2": {"file": "link2.agp", "md5sum": _CONTENT_MD5},
            "link3": {"file": "link3.agp", "md5sum": _CONTENT_MD5},
        }
    }
    maker = ManifestMaker(tmp_path)
    test_files = maker.get_files_checksums()
    assert test_files == expected_content


@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum_warning_subdir(tmp_path: Path, caplog: LogCaptureFixture) -> None:
    """Tests the `ManifestMaker.get_files_checksum()` with a subdir."""
    Path(tmp_path / "sub_dir").mkdir()
    maker = ManifestMaker(tmp_path)
    test_files = maker.get_files_checksums()
    assert caplog.text.strip().endswith("Can't create manifest for subdirectory")
    assert not test_files


@pytest.mark.dependency(depends=["test_init"])
def test_get_files_checksum_warning_empty(tmp_path: Path, caplog: LogCaptureFixture) -> None:
    """Tests the `ManifestMaker.get_files_checksum()` with an empty file, deleted."""
    empty_file = Path(tmp_path / "empty_file.txt")
    empty_file.touch()
    assert empty_file.exists()
    maker = ManifestMaker(tmp_path)
    test_files = maker.get_files_checksums()
    assert "Skip and delete empty file" in caplog.text
    assert not test_files
    assert not empty_file.exists()


@pytest.mark.parametrize(
    "files, expected_content",
    [
        param([], {}, id="No files"),
        param(["genes.gff3"], {"gff3": {"file": "genes.gff3", "md5sum": _CONTENT_MD5}}, id="1 gff3 file"),
    ],
)
@pytest.mark.dependency(depends=["test_init"])
def test_create_manifest(
    tmp_path: Path, assert_files: Callable, files: list[str], expected_content: dict
) -> None:
    """Tests the `ManifestMaker.create_manifest()`."""
    for file_name in files:
        with Path(tmp_path / file_name).open("w") as fh:
            fh.write("CONTENT")

    maker = ManifestMaker(tmp_path)
    maker.create_manifest()
    manifest_file = Path(tmp_path / "manifest.json")
    assert manifest_file.exists()

    expected_manifest_file = Path(tmp_path / "expected.json")
    with expected_manifest_file.open("w") as expected_fh:
        expected_fh.write(json.dumps(expected_content, sort_keys=True, indent=4))
    assert_files(manifest_file, expected_manifest_file)
