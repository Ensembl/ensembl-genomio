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
"""Unit testing of `ensembl.io.genomio.genome_metadata.extend` module.

Typical usage example::
    $ pytest test_extend.py

"""

from pathlib import Path
from typing import List

import pytest

from ensembl.io.genomio.genome_metadata import extend


@pytest.mark.parametrize(
    "gbff_file, output",
    [
        pytest.param("", [], id="No GBFF file"),
        ("input.gbff.gz", ["LR605957", "LR605956"]),
    ],
)
def test_get_gbff_regions(data_dir: Path, gbff_file: str, output: List[str]) -> None:
    """Tests the `extend.get_gbff_regions` class.

    Args:
        data_dir: Module's test data directory fixture.
        gbff_path: GBFF file name.
        output: Expected list of sequence region IDs.
    """
    if gbff_file:
        gbff_path = data_dir / gbff_file
    else:
        gbff_path = None
    result = extend.get_gbff_regions(gbff_path)
    assert result == output
