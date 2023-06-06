#!/usr/bin/env python
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
"""Unit testing of :mod:`ensembl.io.genomio.schemas` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_schemas.py

"""

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager

import pytest

from ensembl.io.genomio.integrity import IntegrityTool


class TestIntegrity:
    """Tests for the integrity module."""

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, tmp_dir: Path):
        """Loads necessary fixtures and values as class attributes."""
        type(self).tmp_dir = tmp_dir

    @pytest.mark.parametrize(
        "manifest_dir, brc_mode, ignore_false_stops, expected",
        [
            (
                [pytest.manifest_dir / "data1/manifest.json"],
                [None, True, False],
                [None, True, False],
                does_not_raise(),
            ),
        ],
    )
    def test_integrity(
        self, manifest_dir: Path, brc_mode: bool, ignore_false_stops: bool, expected: ContextManager
    ) -> None:
        """Tests `integrity:IntegrityTool` method.

        Args:
            brc_mode: BRC specific mode.
            ignore_false_stops: Ignore false stops.

        """
        with expected:
            integrity = IntegrityTool(manifest_dir, brc_mode, ignore_false_stops)
            assert isinstance(integrity, IntegrityTool)

    @pytest.mark.parametrize(
        "manifest_dir, expected",
        [
            (pytest.manifest_dir / "data1/manifest.json", does_not_raise()),
        ],
    )
    def test_get_manifest(self, manifest_dir: Path, expected: ContextManager) -> None:
        """Tests `integrity:IntegrityTool:get_manifest()` method."""
        with expected:
            integrity = IntegrityTool(manifest_dir)
            manifest = integrity.get_manifest()
            assert isinstance(manifest, dict)
