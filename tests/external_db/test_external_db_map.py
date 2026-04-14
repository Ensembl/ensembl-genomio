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
"""Unit testing of `ensembl.io.genomio.external_db.db_map` module."""

from contextlib import nullcontext as no_raise
from pathlib import Path
from typing import ContextManager

import pytest
from pytest import param, raises

from ensembl.io.genomio.external_db.db_map import (
    MapFormatError,
    get_external_db_map,
    DEFAULT_EXTERNAL_DB_MAP,
)


@pytest.mark.parametrize(
    "file_content, expected_output, expected",
    [
        param("", {}, no_raise()),
        param("#Comment", {}, no_raise()),
        param("FOO\tBAR", {"BAR": "FOO"}, no_raise()),
        param("FOO\tBAR\tLOREM", {"BAR": "FOO"}, no_raise()),
        param("FOO", {}, raises(MapFormatError)),
    ],
)
def test_get_external_db_map(
    tmp_path: Path, file_content: str, expected_output: dict, expected: ContextManager
) -> None:
    """Tests the `get_external_db_map` method.

    Args:
        file_content: Test db_map file content.
        expected_output: Db map expected output.
        expected: Context Manager for the expected exception.

    """
    with expected:
        test_file = Path(tmp_path / "test.map")
        with test_file.open("w") as test_fh:
            test_fh.write(file_content)
        output = get_external_db_map(test_file)
        assert output == expected_output


def test_default_map() -> None:
    """Tests the default_map file."""
    with no_raise():
        output = get_external_db_map(DEFAULT_EXTERNAL_DB_MAP)
        assert output
