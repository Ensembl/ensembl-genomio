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
import re
from typing import ContextManager, Optional, Set

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


@pytest.mark.parametrize(
    "seq, pattern, expectation",
    [
        ("", None, [len("")]),
        ("", re.compile("N"), [0]),
        ("CANNAN", None, [len("NANNAN")]),
        ("CANNAN", None, [6]),
        ("CAAAAN", re.compile("N"), [6]),
        ("NAAAAN", re.compile("N"), [1, 6]),
        ("NANNAN", re.compile("N"), [1, 3, 4, 6]),
        ("NANNAN", re.compile("NN"), [4, 6]),
        ("NAAAAN", re.compile("NN"), [6]),
    ],
)
def test_split_seq_by_n(seq: str, pattern: Optional[re.Pattern], expectation: list[int]) -> None:
    """Tests the `chunk.split_seq_by_n` function.

    Args:
        seq: A sequence to split
        pattern: A pattern to split on
        expectation: A list of open chunk ends (like for python list slices)
    """
    assert FastaChunking.split_seq_by_n(seq, pattern) == expectation


@pytest.mark.parametrize(
    "chunk_ends, chunk_size, tolerated_size, expectation",
    [
        (None, 1, None, None),
        ([], 1, None, []),
        ([], 1, 1, []),
        ([6], -1, None, [6]),
        ([6], 6, None, [6]),
        ([6], 5, None, [5, 6]),
        ([4, 6], 0, None, [4, 6]),
        ([4, 6], 1, None, [1, 2, 3, 4, 5, 6]),
        ([4, 6], 2, None, [2, 4, 6]),
        ([4, 6], 2, 4, [4, 6]),
        ([1, 3, 4, 6], 2, 4, [1, 3, 4, 6]),  # do not join back
        ([1, 3, 4, 6], 6, 0, [1, 3, 4, 6]),  # do not join back
    ],
)
def test_split_seq_by_chunk_size(
    chunk_ends: list[int], chunk_size: int, tolerated_size: Optional[int], expectation: list[int]
):
    """Tests the `chunk.split_seq_by_chunk_size` function.

    Args:
        chunk_ends: A list of chunk_ends (open ones, like for python list slices)
        chunk_size: Chunk size
        tolerated_size: A more relaxed value of the chunk size
        expectation: A list of open chunk ends (python slices coordinates)
    """
    assert FastaChunking.split_seq_by_chunk_size(chunk_ends, chunk_size, tolerated_size) == expectation
