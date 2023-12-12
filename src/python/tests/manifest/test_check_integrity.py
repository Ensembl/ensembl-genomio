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
"""Unit testing of `ensembl.io.genomio.manifest.check_integrity` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_check_integrity.py

"""

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager

import pytest

from ensembl.io.genomio.manifest.check_integrity import IntegrityTool, Manifest


class TestIntegrityTool:
    """Tests for the `IntegrityTool` class."""

    tmp_dir: Path
    manifest_dir: Path

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, tmp_dir: Path, manifest_dir: Path):
        """Loads necessary fixtures and values as class attributes."""
        type(self).tmp_dir = tmp_dir
        type(self).manifest_dir = manifest_dir

    @pytest.mark.parametrize(
        "manifest_path, brc_mode, ignore_false_stops, expected",
        [
            (
                "data1/manifest.json",
                [None, True, False],
                [None, True, False],
                does_not_raise(),
            ),
        ],
    )
    def test_check_integrity(
        self, manifest_path: Path, brc_mode: bool, ignore_false_stops: bool, expected: ContextManager
    ) -> None:
        """Tests the `IntegrityTool.check_integrity()` method.

        Args:
            brc_mode: BRC specific mode.
            ignore_false_stops: Ignore false stops.

        """
        manifest_path = self.manifest_dir / manifest_path
        with expected:
            integrity = IntegrityTool(manifest_path, brc_mode, ignore_false_stops)
            assert isinstance(integrity, IntegrityTool)
            assert integrity.brc_mode == brc_mode
            assert integrity.manifest.brc_mode == brc_mode

    @pytest.mark.parametrize(
        "manifest_path, expected",
        [
            ("data1/manifest.json", does_not_raise()),
        ],
    )
    def test_manifest(self, manifest_path: Path, expected: ContextManager) -> None:
        """Tests the `IntegrityTool.manifest` attribute."""
        manifest_dir = self.manifest_dir / manifest_path
        with expected:
            integrity = IntegrityTool(manifest_dir)
            manifest = integrity.manifest
            assert isinstance(manifest, Manifest)
