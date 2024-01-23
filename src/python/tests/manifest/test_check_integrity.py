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


@pytest.mark.parametrize(
    "manifest_file, brc_mode, ignore_false_stops, expected",
    [
        ("manifest.json", True, False, does_not_raise()),
    ],
)
def test_check_integrity(
    data_dir: Path, manifest_file: str, brc_mode: bool, ignore_false_stops: bool, expected: ContextManager
) -> None:
    """Tests the `IntegrityTool.check_integrity()` method.

    Args:
        data_dir: Module's test data directory fixture.
        manifest_file: Manifest file to load.
        brc_mode: BRC specific mode.
        ignore_false_stops: Ignore false stops.
        expected: Context manager for the expected exception, i.e. the test will only pass if that
            exception is raised. Use `contextlib.nullcontext` if no exception is expected.

    """
    with expected:
        integrity = IntegrityTool(data_dir / manifest_file, brc_mode, ignore_false_stops)
        assert integrity.brc_mode == brc_mode
        assert integrity.manifest.brc_mode == brc_mode


@pytest.mark.parametrize(
    "manifest_file, expected",
    [
        ("manifest.json", does_not_raise()),
    ],
)
def test_manifest(data_dir: Path, manifest_file: str, expected: ContextManager) -> None:
    """Tests the `IntegrityTool.manifest` attribute.

    Args:
        data_dir: Module's test data directory fixture.
        manifest_file: Manifest file to load.
        expected: Context manager for the expected exception, i.e. the test will only pass if that
            exception is raised. Use `contextlib.nullcontext` if no exception is expected.

    """
    with expected:
        integrity = IntegrityTool(data_dir / manifest_file)
        assert isinstance(integrity.manifest, Manifest)
