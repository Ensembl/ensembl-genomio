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
"""Unit testing of `ensembl.io.genomio.fasta.chunk` module.

Typical usage example::
    $ pytest test_chunk.py

"""

from contextlib import nullcontext as does_not_raise
import filecmp
from pathlib import Path
from typing import ContextManager, Set

import pytest

import ensembl.io.genomio.fasta.chunk as FastaChunking


@pytest.mark.parametrize(
    "msg, expectation",
    [
        ("A value error test", pytest.raises(ValueError, match=r"^A value error test$")),
    ],
)
def test__on_value_error(msg: str, expectation: ContextManager) -> None:
    """Tests the `chunk._on_value_error` function.

    Args:
        msg: Msg to raise.
        expectation: A context manager with expected exception (`pytest.raises` or nullcontext)
    """

    with expectation:
        FastaChunking._on_value_error(msg)


@pytest.mark.parametrize(
    "chunk_size, chunk_tolerance, expectation",
    [
        (50000, 0, does_not_raise()),
        (50000, 1, does_not_raise()),
        (50001, 1000, does_not_raise()),
        (49999, 0, pytest.raises(ValueError, match=r"wrong.*chunk_size")),
        (49999, -1, pytest.raises(ValueError, match=r"wrong.*chunk_size")),
        (50000, -1, pytest.raises(ValueError, match=r"wrong.*chunk_tolerance")),
    ],
)
def test_check_chunk_size_and_tolerance(
    chunk_size: int,
    chunk_tolerance: int,
    expectation: ContextManager,
) -> None:
    """Tests the `chunk.check_chunk_size_and_tolerance` function.

    Args:
        chunk_size: Chunk size to check
        chunk_tolerance: Chunk tolerance to check
        expectation: A context manager with expected exception (`pytest.raises` or nullcontext)
    """

    with expectation:
        FastaChunking.check_chunk_size_and_tolerance(chunk_size, chunk_tolerance)
