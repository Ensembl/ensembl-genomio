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

from contextlib import contextmanager, nullcontext
from contextlib import nullcontext as does_not_raise
from io import StringIO
import filecmp
from pathlib import Path
import re
from typing import ContextManager, Generator, Optional

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


def test_individual_file_opener(tmp_path: Path) -> None:
    """Tests the `chunk._individual_file_opener` function.

    Args:
        tmp_path: Where temporary files will be created.
    """
    test_dir = Path(tmp_path, "file_opener_test")
    test_dir.mkdir()
    test_path = Path(test_dir, "test.file")

    with FastaChunking._individual_file_opener(str(test_path)) as out:
        print("test", file=out)
        print("test", file=out)
    assert test_path.read_text(encoding="utf-8") == "test\ntest\n"

    with FastaChunking._individual_file_opener(str(test_path)) as out:
        print("test", file=out)
    assert test_path.read_text(encoding="utf-8") == "test\n"

    assert len(list(test_dir.iterdir())) == 1


def test_prepare_out_dir_for_individuals(tmp_path: Path) -> None:
    """Tests the `chunk.prepare_out_dir_for_individuals` function.

    Args:
        tmp_path: Where temporary files will be created.
    """
    test_dir = Path(tmp_path, "prepare_out_dir_test")
    test_file = Path(test_dir, "test.file")

    assert FastaChunking.prepare_out_dir_for_individuals(test_dir, "test.file") == test_file
    assert FastaChunking.prepare_out_dir_for_individuals(test_dir, "test.file") == test_file

    test_file.write_text("test\n", encoding="utf-8")
    assert test_file.read_text(encoding="utf-8") == "test\n"

    assert FastaChunking.prepare_out_dir_for_individuals(test_dir, "test.file") == test_file
    assert len(list(test_dir.iterdir())) == 1


@pytest.mark.parametrize(
    "size, tolerance, expectation",
    [
        (-1, -1, -1),
        (0, -1, 0),
        (10, -1, 10),
        (10, 20, 12),
        (10, 29, 12),
        (10, 30, 13),
        (100, 29, 129),
    ],
)
def test_get_tolerated_size(size: int, tolerance: int, expectation: int) -> None:
    """Tests the `chunk.get_tolerated_size` function.

    Args:
        size: Base size
        tolerance: Percent of allowed deviance as integer.
        expectation: An expexted tolerated size
    """
    assert FastaChunking.get_tolerated_size(size, tolerance) == expectation


@pytest.mark.parametrize(
    "input_fasta_text, chunk_size, chunk_size_tolerated, n_sequece_len, "
    "chunk_sfx, append_offset_to_chunk_name, "
    "expected_chunked_fasta_text, expected_agp_list, expected_individual_files_count, "
    "expected_raised",
    [
        ("", 2, 2, 2, "p", True, "", [], 0, does_not_raise()),
        ("\n", 2, 2, 2, "p", True, "", [], 0, does_not_raise()),
        ("AA\n", 2, 2, 2, "p", True, "", [], 0, does_not_raise()),
        (
            # title, no sequence
            ">a\n",
            2,
            2,
            2,  # Ns
            "p",
            True,
            ">a_p_001_off_0 AGP a 1 0 1 W a_p_001_off_0 1 0 +\n",
            ["a\t1\t0\t1\tW\ta_p_001_off_0\t1\t0\t+"],
            1,
            does_not_raise(),
        ),
        (
            ">c\nAAANNAAA\n",
            3,
            3,
            3,  # Ns
            "p",
            True,
            ">c_p_001_off_0 AGP c 1 3 1 W c_p_001_off_0 1 3 +\n"
            "AAA\n"
            ">c_p_002_off_3 AGP c 4 6 2 W c_p_002_off_3 1 3 +\n"
            "NNA\n"
            ">c_p_003_off_6 AGP c 7 8 3 W c_p_003_off_6 1 2 +\n"
            "AA\n",
            [
                "c\t1\t3\t1\tW\tc_p_001_off_0\t1\t3\t+",
                "c\t4\t6\t2\tW\tc_p_002_off_3\t1\t3\t+",
                "c\t7\t8\t3\tW\tc_p_003_off_6\t1\t2\t+",
            ],
            3,
            does_not_raise(),
        ),
        (
            ">c\nAAANNAAA\n",
            3,
            5,
            2,  # Ns
            "p",
            True,
            ">c_p_001_off_0 AGP c 1 5 1 W c_p_001_off_0 1 5 +\n"
            "AAANN\n"
            ">c_p_002_off_5 AGP c 6 8 2 W c_p_002_off_5 1 3 +\n"
            "AAA\n",
            [
                "c\t1\t5\t1\tW\tc_p_001_off_0\t1\t5\t+",
                "c\t6\t8\t2\tW\tc_p_002_off_5\t1\t3\t+",
            ],
            2,
            does_not_raise(),
        ),
    ],
)
def test_chunk_fasta_stream(
    input_fasta_text: str,
    chunk_size: int,
    chunk_size_tolerated: int,
    n_sequece_len: int,
    chunk_sfx: str,
    append_offset_to_chunk_name: Optional[bool],
    expected_chunked_fasta_text: str,
    expected_agp_list: list[str],
    expected_individual_files_count: int,
    expected_raised: ContextManager,
) -> None:
    """Tests the `chunk.chunk_fasta_stream` function.

    Args:
        tmp_path: Where temporary files will be created.
        ...
    """

    # a workaround for storing individual chunks
    @contextmanager
    def gen_individual_opener(name: str, parts: list[tuple[str, str]]) -> Generator[StringIO]:
        stream = StringIO()
        try:
            yield stream
        finally:
            # Code to release resource, e.g.:
            parts.append((name, stream.getvalue()))
            stream.close()

    parts: list[tuple[str, str]] = []

    def _individual_opener(name: str):
        return gen_individual_opener(name, parts)

    # assert joined case
    with StringIO(input_fasta_text) as input_fasta:
        with StringIO() as output_fasta:
            parts.clear()
            with expected_raised:
                agp_list = FastaChunking.chunk_fasta_stream(
                    input_fasta,
                    chunk_size,
                    chunk_size_tolerated,
                    output_fasta,
                    None,
                    n_sequece_len,
                    chunk_sfx,
                    append_offset_to_chunk_name,
                    open_individual=_individual_opener,
                )
            assert output_fasta.getvalue() == expected_chunked_fasta_text
            assert agp_list == expected_agp_list
            assert len(parts) == 0

    with StringIO(input_fasta_text) as input_fasta:
        with nullcontext() as no_output_fasta:
            parts.clear()
            with expected_raised:
                agp_list = FastaChunking.chunk_fasta_stream(
                    input_fasta,
                    chunk_size,
                    chunk_size_tolerated,
                    no_output_fasta,
                    "/path/prefix",
                    n_sequece_len,
                    chunk_sfx,
                    append_offset_to_chunk_name,
                    open_individual=_individual_opener,
                )
            assert "".join(map(lambda p: p[1], parts)) == expected_chunked_fasta_text
            assert agp_list == expected_agp_list
            assert len(parts) == expected_individual_files_count


# assert individual files case

# assert FastaChunking.chunk_fasta_stream(...) == expected_agp

# and for inividual files output
#    check joined files == joined_output
#    check for file count
