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
"""Unit testing of `ensembl.io.genomio.manifest.check_integrity` module."""

from pathlib import Path

import pytest

from ensembl.io.genomio.manifest.check_integrity import IntegrityTool, ManifestStats


@pytest.mark.parametrize(
    ("manifest_file", "ignore_false_stops"),
    [
        ("manifest.json", False),
        ("manifest.json", True),
    ],
)
def test_check_integrity(
    data_dir: Path,
    manifest_file: str,
    ignore_false_stops: bool,
) -> None:
    """Test the `IntegrityTool.check_integrity()` method.

    Args:
        data_dir: Module's test data directory fixture.
        manifest_file: Manifest file to load.
        ignore_false_stops: Ignore false stops.

    """
    integrity = IntegrityTool(data_dir / manifest_file, ignore_false_stops)
    assert integrity.ignore_final_stops == ignore_false_stops


@pytest.mark.parametrize(
    "manifest_file",
    ["manifest.json"],
)
def test_manifest(data_dir: Path, manifest_file: str) -> None:
    """Test the `IntegrityTool.manifest` attribute.

    Args:
        data_dir: Module's test data directory fixture.
        manifest_file: Manifest file to load.

    """
    integrity = IntegrityTool(data_dir / manifest_file)
    assert isinstance(integrity.manifest, ManifestStats)
