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
"""Unit testing of `ensembl.io.genomio.seq_region.dump` module."""

from contextlib import nullcontext as no_raise
from pathlib import Path
from typing import ContextManager

import pytest
from pytest import param, raises

from ensembl.io.genomio.seq_region.dump import fetch_coord_systems

def test_fetch_coord_systems(
) -> None:
    """Tests the `fetch_coord_system` method."""
    assert True
