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

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager

import pytest
from pytest import param

from ensembl.io.genomio.manifest.generate import ManifestMaker


@pytest.mark.parametrize(
    "files, expected_manifest_files",
    [
        param({}, [], id="No files"),
        param({"foobar.json", "manifest.json"}, [], id="No recognizable file"),
    ],
)
def test_get_files_checksum(tmp_path: Path, files: list[str], expected_manifest_files: dict) -> None:
    """Tests the `ManifestMaker.get_files_checksum()` method.

    Args:
        tmp_path: Test tmp dir.
        files: Name of files to create for the test.
        expected_manifest_files: Dict of manifest files loaded.

    """
    maker = ManifestMaker(tmp_path)
    for file in files:
        with Path(tmp_path / file).open("w") as fh:
            fh.write("CONTENT")
    test_files = maker.get_files_checksums()
    assert set(test_files) == set(expected_manifest_files)
