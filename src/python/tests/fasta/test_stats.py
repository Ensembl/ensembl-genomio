# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Unit testing of `ensembl.io.genomio.fasta.stats` module."""

import gzip
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from pytest import param

from ensembl.io.genomio.fasta import stats


def write_text(path: Path, text: str) -> Path:
    """Write text to a test file."""
    path.write_text(text, encoding="utf-8")
    return path


def write_gzip_text(path: Path, text: str) -> Path:
    """Write gzip-compressed text to a test file."""
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    return path


@pytest.mark.parametrize(
    "contents, expected",
    [
        param(">seq1\nACGT\nAC\n>seq2 description\nNNN\n", "6 9 2\n", id="Multiline records"),
        param(">empty\n>seq\nAC\n", "2 2 2\n", id="Empty record"),
        param("", "0 0 0\n", id="Empty file"),
        param("ACGT\n", "0 0 0\n", id="Sequence without record header"),
    ],
)
def test_compute_fasta_stats(tmp_path: Path, contents: str, expected: str) -> None:
    """
    Tests the `stats.compute_fasta_stats()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        contents: FASTA file content.
        expected: Expected stats file content.
    """
    fasta_file = write_text(tmp_path / "in.fa", contents)
    output_file = tmp_path / "stats.txt"

    assert stats.compute_fasta_stats(fasta_file=fasta_file, output_file=output_file) is None
    assert output_file.read_text(encoding="utf-8") == expected


def test_compute_fasta_stats_default_output_file(tmp_path: Path) -> None:
    """
    Tests the `stats.compute_fasta_stats()` function with the default output file.

    Args:
        tmp_path: Test's unique temporary directory fixture.
    """
    fasta_file = write_text(tmp_path / "in.fa", ">seq\nACGT\n")
    output_file = tmp_path / "in.stats.txt"

    stats.compute_fasta_stats(fasta_file=fasta_file, output_file=None)

    assert output_file.read_text(encoding="utf-8") == "4 4 1\n"


def test_compute_fasta_stats_reads_gzipped_fasta(tmp_path: Path) -> None:
    """
    Tests the `stats.compute_fasta_stats()` function with gzipped input.

    Args:
        tmp_path: Test's unique temporary directory fixture.
    """
    fasta_file = write_gzip_text(tmp_path / "in.fa.gz", ">seq\nAAAA\n")
    output_file = tmp_path / "stats.txt"

    stats.compute_fasta_stats(fasta_file=fasta_file, output_file=output_file)

    assert output_file.read_text(encoding="utf-8") == "4 4 1\n"


def test_write_fasta_stats(tmp_path: Path) -> None:
    """
    Tests the `stats.write_fasta_stats()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
    """
    output_file = tmp_path / "stats.txt"

    stats.write_fasta_stats(stats.FastaStats(longest=4, total=10, n_seqs=3), output_file)

    assert output_file.read_text(encoding="utf-8") == "4 10 3\n"


@pytest.mark.parametrize(
    "argv, expected_output",
    [
        param(["--fasta"], None, id="Default output"),
        param(["--fasta", "--output"], "stats.txt", id="Explicit output"),
    ],
)
def test_parse_args(tmp_path: Path, argv: list[str], expected_output: str | None) -> None:
    """
    Tests the `stats.parse_args()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        argv: Argument names to pass to the parser.
        expected_output: Expected output argument value.
    """
    fasta_file = write_text(tmp_path / "in.fa", ">seq\nA\n")
    args_list = [argv[0], str(fasta_file)]
    output_file = None
    if expected_output is not None:
        output_file = tmp_path / expected_output
        args_list.extend([argv[1], str(output_file)])

    args = stats.parse_args(args_list)

    assert args.fasta == fasta_file
    if expected_output is None:
        assert not hasattr(args, "output")
    else:
        assert args.output == output_file


@pytest.mark.parametrize(
    "argv, expected_output",
    [
        param(["--fasta"], "in.stats.txt", id="Default output"),
        param(["--fasta", "--output"], "stats.txt", id="Explicit output"),
    ],
)
def test_main(tmp_path: Path, argv: list[str], expected_output: str) -> None:
    """
    Tests the `stats.main()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        argv: Argument names to pass to the parser.
        expected_output: Expected stats output file name.
    """
    fasta_file = write_text(tmp_path / "in.fa", ">seq\nACGT\n")
    args_list = [argv[0], str(fasta_file)]
    if len(argv) > 1:
        args_list.extend([argv[1], str(tmp_path / expected_output)])

    assert stats.main(args_list) is None
    assert (tmp_path / expected_output).read_text(encoding="utf-8") == "4 4 1\n"


@patch.object(stats, "compute_fasta_stats")
def test_main_calls_compute_fasta_stats(mock_compute_fasta_stats: Mock, tmp_path: Path) -> None:
    """
    Tests that `stats.main()` calls `stats.compute_fasta_stats()` with parsed arguments.

    Args:
        mock_compute_fasta_stats: Mock object for the `stats.compute_fasta_stats()` function.
        tmp_path: Test's unique temporary directory fixture.
    """
    fasta_file = write_text(tmp_path / "in.fa", ">seq\nA\n")
    output_file = tmp_path / "stats.txt"

    stats.main(["--fasta", str(fasta_file), "--output", str(output_file)])

    mock_compute_fasta_stats.assert_called_once_with(fasta_file=fasta_file, output_file=output_file)


@patch.object(stats, "compute_fasta_stats")
def test_main_raise_exception(mock_compute_fasta_stats: Mock, tmp_path: Path) -> None:
    """
    Tests the `stats.main()` function when FASTA stats computation fails.

    Args:
        mock_compute_fasta_stats: Mock object for the `stats.compute_fasta_stats()` function.
        tmp_path: Test's unique temporary directory fixture.
    """
    fasta_file = write_text(tmp_path / "in.fa", ">seq\nA\n")
    mock_compute_fasta_stats.side_effect = RuntimeError("Mocked exception")

    with pytest.raises(RuntimeError, match=r"Mocked exception"):
        stats.main(["--fasta", str(fasta_file)])
