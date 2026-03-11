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
"""Unit testing of `ensembl.io.genomio.fasta.split` module."""


from contextlib import nullcontext as does_not_raise
import filecmp
from io import TextIOWrapper
from pathlib import Path
import re
from typing import Any, Callable, ContextManager
from unittest.mock import Mock, patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from deepdiff import DeepDiff
import pytest
from pytest import MonkeyPatch, param

from ensembl.io.genomio.fasta import split


def force_open_failure_for_suffix(suffix: str) -> Callable:
    """Monkeypatches `open` to raise an OSError when trying to open a file with the given suffix."""
    real_open = open

    def _patched_open(self: Path, *args: Any, **kwargs: Any) -> TextIOWrapper | Exception:
        if self.suffix == suffix:
            raise OSError(f"Simulated open failure for files ending with '{suffix}'")
        return real_open(self, *args, **kwargs)

    return _patched_open


@pytest.mark.parametrize(
    "name,expected",
    [("in.fa", "in"), ("in.fa.gz", "in"), ("in", "in")],
)
def test_get_fasta_basename(tmp_path: Path, name: str, expected: str) -> None:
    """
    Tests the `split._get_fasta_basename()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        name: File name.
        expected: Expected value returned by the function.

    """
    tmp_file = tmp_path / name
    tmp_file.touch()
    assert split._get_fasta_basename(tmp_file) == expected


@pytest.mark.parametrize(
    "tree, expectation",
    [
        (Path("1/1/in.1.fa"), does_not_raise()),
        (Path("1/1/unexpected.gff3"), pytest.raises(RuntimeError, match=r"Unexpected file identified")),
    ],
)
def test_check_contents_deletable(tmp_path: Path, tree: Path, expectation: ContextManager) -> None:
    """
    Tests the `split._check_contents_deletable()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        tree: Path to the file tree to create under the temporary directory.
        expectation: Context manager for the expected exception. Use `~contextlib.nullcontext` with
            the expected output if no exception is expected.

    """
    out_dir = tmp_path / "out"
    out_tree = out_dir / tree
    out_tree.parent.mkdir(parents=True, exist_ok=True)
    out_tree.touch()
    output_file_re = re.compile(r"^in\.[1-9]\d*(\.[1-9]\d*)?\.fa$")
    with expectation:
        split._check_contents_deletable(out_dir, output_file_re)
    assert out_tree.exists()


@pytest.mark.parametrize(
    "fasta_file, tree, expected",
    [
        param(Path("in.fa"), [], [], id="No output dir"),
        param(
            Path("in.fa"),
            [Path("1/in.1.fa"), Path("2/in.2.fa"), Path("keep")],
            [False, False, True],
            id="Default case without AGP file",
        ),
        param(Path("in.fa"), [Path("in.agp")], [False], id="Delete AGP file"),
    ],
)
def test_clean_previous_output(tmp_path: Path, fasta_file: Path, tree: list[Path], expected: str) -> None:
    """
    Tests the `split._clean_previous_output()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        fasta_file: Path to the input FASTA file.
        tree: List of output files to create for the test.
        expected: List of whether the corresponding file in the tree is expected to exist after cleaning.

    """
    out_dir = tmp_path / "out"
    for rel_path in tree:
        out_path = out_dir / rel_path
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.touch()
    split._clean_previous_output(fasta_file, out_dir)
    if not tree:
        assert not out_dir.exists()
    else:
        assert out_dir.exists()
        for i, rel_path in enumerate(tree):
            assert (out_dir / rel_path).exists() == expected[i]
            # Expect folders inside out_dir to be deleted as well
            if rel_path.parent != Path("."):
                assert (out_dir / rel_path.parent).exists() == expected[i]


class TestOutputWriter:
    """Tests `split.OutputWriter` class."""

    @pytest.mark.parametrize(
        "fasta_file, write_agp, unique_file_names, max_files, max_dirs, expected_out_path, expected_agp_name",
        [
            param(Path("in.fa"), False, False, None, None, "1/in.1.fa", "", id="Default args"),
            param(Path("in.fa"), True, False, None, None, "1/in.1.fa", "in.agp", id="AGP enabled"),
            param(Path("in.fa"), False, True, None, None, "1/in.0.1.fa", "", id="Unique file names"),
            param(Path("in.fa"), False, False, 2, 2, "1/in.1.fa", "", id="Set max elements per dir"),
        ],
    )
    def test_init(
        self,
        tmp_path: Path,
        fasta_file: Path,
        write_agp: bool,
        unique_file_names: bool,
        max_files: int | None,
        max_dirs: int | None,
        expected_out_path: str,
        expected_agp_name: str,
    ) -> None:
        """
        Tests the `__init__()` method of the `split.OutputWriter` class.

        Args:
            tmp_path: Test's unique temporary directory fixture.
            fasta_file: Input raw or compressed FASTA file containing sequences to split.
            write_agp: Write an AGP v2.0 file describing how each input sequence maps to output chunks.
            unique_file_names: Include folder index in output FASTA filenames to make them unique.
            max_files_per_directory: Maximum number of FASTA files per directory.
            max_dirs_per_directory: Maximum number of subdirectories per directory level.
            expected_out_path: Expected relative path to the output FASTA file.
            expected_agp_name: Expected name of the AGP file in the output directory.

        """
        out_dir = tmp_path / "out"
        writer = split.OutputWriter(
            fasta_file=fasta_file,
            out_dir=out_dir,
            write_agp=write_agp,
            unique_file_names=unique_file_names,
            max_files_per_directory=max_files,
            max_dirs_per_directory=max_dirs,
        )
        writer.close()
        assert writer
        assert writer.basename == fasta_file.stem
        assert writer.file_count == 1
        assert (out_dir / expected_out_path).exists()
        if write_agp:
            assert writer.agp_file == out_dir / expected_agp_name and writer.agp_file.exists()
            with writer.agp_file.open("r") as agp_fh:
                header = agp_fh.readline().strip()
                assert header == "# AGP-version 2.0"
        else:
            assert writer.agp_file is None

    @pytest.mark.parametrize(
        "write_agp, suffix, exc_msg",
        [
            param(False, ".fa", "Failed to open output file", id="FASTA file"),
            param(True, ".agp", "Failed to open AGP file", id="AGP file"),
        ],
    )
    def test_create_file_exception(
        self, monkeypatch: MonkeyPatch, tmp_path: Path, write_agp: bool, suffix: str, exc_msg: str
    ) -> None:
        """
        Tests the `_create_output_file()` and `_create_agp_file()` methods of the `split.OutputWriter`
        class when an OSError is raised.

        Args:
            monkeypatch: Pytest fixture to patch methods.
            tmp_path: Test's unique temporary directory fixture.
            write_agp: Write an AGP v2.0 file describing how each input sequence maps to output chunks.
            suffix: Suffix of the file for which to simulate the open failure.
            exc_msg: Expected message in the raised RuntimeError.

        """
        monkeypatch.setattr("pathlib.Path.open", force_open_failure_for_suffix(suffix))
        with pytest.raises(RuntimeError, match=rf"{exc_msg}"):
            writer = split.OutputWriter(
                fasta_file=Path("in.fa"), out_dir=tmp_path, write_agp=write_agp, unique_file_names=False
            )
            writer.close()

    def test_open_new_file(self, tmp_path: Path) -> None:
        """
        Tests the `open_new_file()` method of the `split.OutputWriter` class.

        Args:
            tmp_path: Test's unique temporary directory fixture.

        """
        out_dir = tmp_path / "out"
        writer = split.OutputWriter(
            fasta_file=Path("in.fa"), out_dir=out_dir, write_agp=False, unique_file_names=False
        )
        writer.open_new_file()
        writer.close()
        assert writer.file_count == 2
        assert (out_dir / "1" / "in.1.fa").exists()
        assert (out_dir / "1" / "in.2.fa").exists()

    @pytest.mark.parametrize(
        "write_agp, agp_obj_id, agp_start, agp_end, agp_part_nr, expectation",
        [
            param(False, None, None, None, None, does_not_raise(), id="Default args"),
            param(True, "seq1", 1, 4, 1, does_not_raise(), id="Write AGP"),
            param(
                True,
                None,
                1,
                4,
                1,
                pytest.raises(ValueError, match=r"All AGP arguments must be provided if writing AGP entries"),
                id="Missing AGP args",
            ),
        ],
    )
    def test_write_record(
        self,
        tmp_path: Path,
        write_agp: bool,
        agp_obj_id: str | None,
        agp_start: int | None,
        agp_end: int | None,
        agp_part_nr: int | None,
        expectation: ContextManager,
    ) -> None:
        """
        Tests the `write_record()` method of the `split.OutputWriter` class.

        Args:
            tmp_path: Test's unique temporary directory fixture.
            write_agp: Write an AGP v2.0 file describing how each input sequence maps to output chunks.
            agp_obj_id: Original (unchunked) input sequence ID to write into the AGP ``object`` column.
            agp_start: Start coordinate on the AGP object (1-based, inclusive).
            agp_end: End coordinate on the AGP object (1-based, inclusive).
            agp_part_nr: Component part number for this object (starts at 1 per object).
            expectation: Context manager for the expected exception. Use `~contextlib.nullcontext` with
                the expected output if no exception is expected.

        """
        out_dir = tmp_path / "out"
        writer = split.OutputWriter(
            fasta_file=Path("in.fa"), out_dir=out_dir, write_agp=write_agp, unique_file_names=False
        )
        assert writer.record_count == 0
        assert writer.file_len == 0
        in_record = SeqRecord(Seq("ACGT"), id="seq1", description="test sequence")
        with expectation:
            writer.write_record(
                record=in_record,
                agp_object_id=agp_obj_id,
                agp_start=agp_start,
                agp_end=agp_end,
                agp_part_nr=agp_part_nr,
            )
            writer.close()
            out_fasta_file = out_dir / "1" / "in.1.fa"
            assert writer.record_count == 1
            assert writer.file_len == 4
            out_record = SeqIO.read(out_fasta_file, "fasta")
            assert out_record.id == in_record.id
            assert out_record.seq == in_record.seq
            assert out_record.description == f"{in_record.id} {in_record.description}"
            if write_agp:
                with writer.agp_file.open("r") as agp_fh:  # type: ignorep[union-attr]
                    # Skip header line
                    agp_fh.readline()
                    agp_line = agp_fh.readline().strip()
                    expected_line = (
                        f"{agp_obj_id}\t{agp_start}\t{agp_end}\t{agp_part_nr}\tW\t"
                        f"{out_record.id}\t1\t{len(out_record.seq)}\t+"
                    )
                    assert agp_line == expected_line


@pytest.mark.parametrize(
    "extra_args, expected",
    [
        ({}, "default"),
        ({"max_seqs_per_file": 1}, "1_seq"),
        ({"max_seq_length_per_file": 6}, "1_seq"),
        ({"max_seq_length_per_file": 6, "force_max_seq_length": True}, "6bp_force"),
        (
            {"max_seq_length_per_file": 6, "force_max_seq_length": True, "min_chunk_length": 4},
            "6bp_force_min_chunk",
        ),
        (
            {"max_seq_length_per_file": 4, "force_max_seq_length": True, "min_chunk_length": 4},
            "4bp_force_min_chunk",
        ),
    ],
)
def test_split_fasta(tmp_path: Path, data_dir: Path, extra_args: dict[str, Any], expected: str) -> None:
    """
    Tests the `split.split_fasta()` function.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        data_dir: Module's test data directory fixture.
        extra_args:
        expected:

    """
    in_fasta = data_dir / "input.fa"
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)
    split.split_fasta(in_fasta, out_dir, **extra_args)
    report = filecmp.dircmp(data_dir / expected, out_dir)
    report.subdirs["."] = report
    for diff_report in report.subdirs.values():
        assert diff_report.left_only == []
        assert diff_report.right_only == []
        assert diff_report.diff_files == []


def test_split_fasta_empty_file(tmp_path: Path) -> None:
    """
    Tests the `split.split_fasta()` function when an empty input file is provided.

    Args:
        tmp_path: Test's unique temporary directory fixture.

    """
    in_fasta = tmp_path / "empty.fa"
    in_fasta.touch()
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)
    split.split_fasta(in_fasta, out_dir)
    assert list(out_dir.iterdir()) == []


def test_split_fasta_rm_existing_files(tmp_path: Path, data_dir: Path) -> None:
    """
    Tests the `split.split_fasta()` function when there are files from a previous run and we want to
    delete them.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        data_dir: Module's test data directory fixture.

    """
    in_fasta = data_dir / "input.fa"
    out_dir = tmp_path / "out"
    # Create a subtree that resembles what would have been produced in a previous run
    (out_dir / "2").mkdir(parents=True, exist_ok=True)
    (out_dir / "2" / "input.2.fa").touch()
    split.split_fasta(in_fasta, out_dir, delete_existing_files=True)
    report = filecmp.dircmp(data_dir / "default", out_dir)
    report.subdirs["."] = report
    for diff_report in report.subdirs.values():
        assert diff_report.left_only == []
        assert diff_report.right_only == []
        assert diff_report.diff_files == []


@pytest.mark.parametrize(
    "arg_list, expectation",
    [
        param(
            ["--fasta-file", __file__],
            does_not_raise(
                {
                    "fasta_file": __file__,
                    "write_agp": False,
                    "delete_existing_files": False,
                    "unique_file_names": False,
                    "force_max_seq_length": False,
                    "log_level": "WARNING",
                },
            ),
            id="Default args",
        ),
        param(
            [
                "--fasta-file",
                __file__,
                "--out-dir",
                str(Path(__file__).parent),
                "--write-agp",
                "--delete-existing-files",
                "--unique-file-names",
                "--force-max-seq-length",
                "--max-seqs-per-file",
                "1",
                "--max-seq-length-per-file",
                "1",
                "--min-chunk-length",
                "1",
                "--max-files-per-directory",
                "1",
                "--max-dirs-per-directory",
                "1",
            ],
            does_not_raise(
                {
                    "fasta_file": __file__,
                    "out_dir": str(Path(__file__).parent),
                    "write_agp": True,
                    "delete_existing_files": True,
                    "unique_file_names": True,
                    "force_max_seq_length": True,
                    "max_seqs_per_file": 1,
                    "max_seq_length_per_file": 1,
                    "min_chunk_length": 1,
                    "max_files_per_directory": 1,
                    "max_dirs_per_directory": 1,
                    "log_level": "WARNING",
                },
            ),
            id="New arg values",
        ),
        param(
            ["--fasta-file", __file__, "--min-chunk-length", "2"],
            pytest.raises(ValueError, match=r"--min-chunk-length requires --max-seq-length-per-file"),
            id="min_chunk_length without max_seq_length_per_file",
        ),
    ],
)
def test_parse_args(arg_list: list[str], expectation: ContextManager) -> None:
    """
    Tests the `split.parse_args()` function.

    Args:
        arg_list: List of command line arguments to parse.
        expectation: Context manager for the expected exception. Use `~contextlib.nullcontext` with
            the expected output if no exception is expected.
    """
    with expectation as exp:
        args = split.parse_args(arg_list)
        # DeepDiff is not able to compare two objects of Path type - need to convert them to string
        setattr(args, "fasta_file", str(args.fasta_file))
        if hasattr(args, "out_dir"):
            setattr(args, "out_dir", str(args.out_dir))
        assert not DeepDiff(vars(args), exp)


@patch("ensembl.io.genomio.fasta.split.split_fasta")
def test_main(mock_split_fasta: Mock, tmp_path: Path) -> None:
    """
    Tests the `split.main()` function (entry point).

    Args:
        mock_split_fasta: Mock object for the `split.split_fasta()` function.
        tmp_path: Temporary directory provided by pytest.
    """
    fasta_path = tmp_path / "in.fa"
    fasta_path.touch()
    split.main(["--fasta-file", str(fasta_path)])
    # Check that we have called the mocked function once with the expected parameters
    mock_split_fasta.assert_called_once_with(
        fasta_file=fasta_path,
        out_dir=None,
        write_agp=False,
        delete_existing_files=False,
        unique_file_names=False,
        force_max_seq_length=False,
        max_seqs_per_file=None,
        max_seq_length_per_file=None,
        min_chunk_length=None,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )


@patch("ensembl.io.genomio.fasta.split.split_fasta")
def test_main_raise_exception(mock_split_fasta: Mock, tmp_path: Path) -> None:
    """
    Tests the `split.main()` function (entry point).

    Args:
        mock_split_fasta: Mock object for the `split.split_fasta()` function.
        tmp_path: Temporary directory provided by pytest.
    """
    fasta_path = tmp_path / "in.fa"
    fasta_path.touch()
    mock_split_fasta.side_effect = RuntimeError("Mocked exception")
    with pytest.raises(RuntimeError, match=r"Mocked exception"):
        split.main(["--fasta-file", str(fasta_path)])
