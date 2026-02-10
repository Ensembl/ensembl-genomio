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
import importlib
from pathlib import Path
import re

import pytest
from Bio import SeqIO

from .utils import read_fasta


@pytest.fixture(scope="session")
def fasta_split():
    return importlib.import_module("ensembl.io.genomio.fasta.split")


def list_output_fastas(out_dir: Path) -> list[Path]:
    return sorted(p for p in out_dir.rglob("*.fa") if p.is_file())


def read_agp_lines(agp_path: Path) -> list[str]:
    return agp_path.read_text(encoding="utf-8").splitlines()


@pytest.mark.parametrize(
    "name,expected",
    [("in.fa", "in"), ("in.fa.gz", "in"), ("in", "in")],
)
def test_get_fasta_basename(fasta_split, tmp_path, name, expected):
    p = tmp_path / name
    p.write_text(">x\nA\n", encoding="utf-8")
    assert fasta_split._get_fasta_basename(p) == expected


def test_clean_previous_output_deletes_numeric_top_level_dirs(fasta_split, tmp_path, write_fasta):
    in_fa = write_fasta("in.fa", [("seq1", "ACGT", None)])
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    # Create numeric dirs with valid outputs
    d1 = out_dir / "1"
    d2 = out_dir / "2"
    d1.mkdir()
    d2.mkdir()
    (d1 / "in.1.fa").write_text(">x\nA\n", encoding="utf-8")
    (d2 / "in.2.1.fa").write_text(">y\nC\n", encoding="utf-8")  # accept both naming modes

    # Non-numeric dir should be untouched, as should numeric folders with preceding zeroes or named '0'
    other = out_dir / "misc"
    other.mkdir()
    (other / "in.3.fa").write_text(">other\nT\n", encoding="utf-8")
    other_numeric = out_dir / "01"
    other_numeric.mkdir()
    (other_numeric / "in.4.fa").write_text(">other_numeric\nG\n", encoding="utf-8")
    other_zero = out_dir / "0"
    other_zero.mkdir()
    (other_zero / "in.5.fa").write_text(">other_zero\nA\n", encoding="utf-8")

    fasta_split.clean_previous_output(in_fa, out_dir)

    assert not d1.exists()
    assert not d2.exists()
    assert other.exists()
    assert other_numeric.exists()
    assert other_zero.exists()
    assert (other / "in.3.fa").exists()
    assert (other_numeric / "in.4.fa").exists()
    assert (other_zero / "in.5.fa").exists()


def test_clean_previous_output_aborts_on_unexpected_file_and_deletes_nothing(
    fasta_split, tmp_path, write_fasta
):
    in_fa = write_fasta("in.fa", [("seq1", "ACGT", None)])
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    d1 = out_dir / "1"
    d2 = out_dir / "2"
    d1.mkdir()
    d2.mkdir()
    (d1 / "in.1.fa").write_text(">x\nA\n", encoding="utf-8")
    (d2 / "unexpected.fa").write_text(">x\nU\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match=r"Unexpected file identified.*out/2/unexpected.fa"):
        fasta_split.clean_previous_output(in_fa, out_dir)

    # Nothing deleted because validation failed before deletion loop
    assert d1.exists()
    assert d2.exists()
    assert (d1 / "in.1.fa").exists()
    assert (d2 / "unexpected.fa").exists()


def test_outputwriter_basename_and_first_file_created(fasta_split, tmp_path):
    in_fa = tmp_path / "in.fa.gz"
    in_fa.write_text(">x\nACGT\n", encoding="utf-8")

    out = tmp_path / "out"
    w = fasta_split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        unique_file_names=False,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )
    try:
        assert w.basename == "in"
        fastas = list_output_fastas(out)
        # open_new_file called in __init__ -> one file exists
        assert len(fastas) == 1
        assert fastas[0].name == "in.1.fa"
    finally:
        w.close()


@pytest.mark.parametrize(
    "max_dirs_per_directory,dir_index,expected",
    [
        (None, 1, ["1"]),
        (None, 99, ["1"]),
        (10, 1, ["1"]),
        (10, 2, ["2"]),
        (10, 10, ["10"]),
        (10, 11, ["1", "1"]),
        (10, 12, ["1", "2"]),
        (2, 1, ["1"]),
        (2, 2, ["2"]),
        (2, 3, ["1", "1"]),
        (2, 4, ["1", "2"]),
    ],
)
def test_get_subdir_path_math(
    fasta_split, write_fasta, tmp_path, max_dirs_per_directory, dir_index, expected
):
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = fasta_split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        unique_file_names=False,
        max_files_per_directory=1,
        max_dirs_per_directory=max_dirs_per_directory,
    )
    try:
        sub = writer._get_subdir_path(dir_index)
        # compare only relative parts under out
        relative = sub.relative_to(out)
        assert list(relative.parts) == expected
    finally:
        writer.close()


def test_file_and_dir_index_with_max_files(fasta_split, write_fasta, tmp_path):
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = fasta_split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        unique_file_names=False,
        max_files_per_directory=3,
        max_dirs_per_directory=10,
    )
    try:
        # file_count starts at 1 after open_new_file in __init__
        assert writer._get_file_and_dir_index() == (1, 1)

        writer._get_path_for_next_file()  # file_count -> 2
        assert writer._get_file_and_dir_index() == (2, 1)

        writer._get_path_for_next_file()  # 3
        assert writer._get_file_and_dir_index() == (3, 1)

        writer._get_path_for_next_file()  # 4 -> rollover to dir_index=2, file_index=1
        assert writer._get_file_and_dir_index() == (1, 2)
    finally:
        writer.close()


def test_unique_file_names_include_dir_index(fasta_split, write_fasta, tmp_path):
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = fasta_split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        unique_file_names=True,
        max_files_per_directory=1,
        max_dirs_per_directory=10,
    )
    try:
        # First file created in __init__
        f1 = list_output_fastas(out)[0]
        assert re.match(r"in\.1\.1\.fa$", f1.name)

        writer.open_new_file()
        f2 = list_output_fastas(out)[-1]
        assert re.match(r"in\.2\.1\.fa$", f2.name)
    finally:
        writer.close()


def test_add_agp_entry_when_agp_disabled(fasta_split, write_fasta, tmp_path):
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = fasta_split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        unique_file_names=False,
    )
    try:
        # should not raise; _agp_fh is None
        writer.add_agp_entry("obj", 1, 1, 1, "part", 1)
        assert not (out / "in.agp").exists()
    finally:
        writer.close()


def test_split_fasta_empty_input_no_outputs(fasta_split, tmp_path):
    in_fa = tmp_path / "empty.fa"
    in_fa.write_bytes(b"")
    out = tmp_path / "out"

    fasta_split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
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
    assert not out.exists() or len(list(out.rglob("*"))) == 0


def test_split_by_max_seqs_per_file_rollover(fasta_split, tmp_path, write_fasta):
    in_fa = write_fasta(
        "in.fa",
        [("a", "AA", None), ("b", "CC", None), ("c", "GG", None), ("d", "TT", None), ("e", "AT", None)],
    )
    out = tmp_path / "out"

    fasta_split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        delete_existing_files=False,
        unique_file_names=False,
        force_max_seq_length=False,
        max_seqs_per_file=2,
        max_seq_length_per_file=None,
        min_chunk_length=None,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )

    fastas = list_output_fastas(out)
    assert len(fastas) == 3

    counts = [len(read_fasta(fasta)) for fasta in fastas]
    assert counts == [2, 2, 1]


def test_split_by_max_seq_length_per_file_rollover(fasta_split, tmp_path, write_fasta):
    in_fa = write_fasta("in.fa", [("a", "AAAA", None), ("b", "TTTT", None)])
    out = tmp_path / "out"

    fasta_split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        delete_existing_files=False,
        unique_file_names=False,
        force_max_seq_length=False,
        max_seqs_per_file=None,
        max_seq_length_per_file=4,
        min_chunk_length=None,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )

    fastas = list_output_fastas(out)
    assert len(fastas) == 2


def test_force_chunking_splits_long_record(fasta_split, tmp_path, write_fasta):
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", "desc")])
    out = tmp_path / "out"

    fasta_split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        delete_existing_files=False,
        unique_file_names=False,
        force_max_seq_length=True,
        max_seqs_per_file=None,
        max_seq_length_per_file=4,
        min_chunk_length=None,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )

    fastas = list_output_fastas(out)
    ids = []
    seqs = []
    for fasta in fastas:
        for record_id, seq in read_fasta(fasta).items():
            ids.append(record_id)
            seqs.append(seq)
    assert ids == ["X_chunk_start_0", "X_chunk_start_4", "X_chunk_start_8"]
    assert seqs == ["ATCG", "GATT", "AC"]


def test_force_chunking_merges_small_remainder(fasta_split, tmp_path, write_fasta, caplog):
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", None)])
    out = tmp_path / "out"

    fasta_split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        delete_existing_files=False,
        unique_file_names=False,
        force_max_seq_length=True,
        max_seqs_per_file=None,
        max_seq_length_per_file=4,
        min_chunk_length=3,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )

    assert any("merging with previous chunk" in r.message for r in caplog.records)

    fastas = list_output_fastas(out)
    ids = []
    seqs = []
    for fasta in fastas:
        for record_id, seq in read_fasta(fasta).items():
            ids.append(record_id)
            seqs.append(seq)

    assert ids == ["X_chunk_start_0", "X_chunk_start_4"]
    assert seqs == ["ATCG", "GATTAC"]


def test_overlong_record_without_chunking_warns_and_written_whole(fasta_split, tmp_path, write_fasta, caplog):
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", None)])
    out = tmp_path / "out"

    fasta_split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        delete_existing_files=False,
        unique_file_names=False,
        force_max_seq_length=False,
        max_seqs_per_file=None,
        max_seq_length_per_file=4,
        min_chunk_length=None,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )

    assert any("chunking not enabled" in r.message for r in caplog.records)

    fastas = list_output_fastas(out)
    assert len(fastas) == 1

    record_id, seq = next(iter(read_fasta(fastas[0]).items()))
    assert record_id == "X"
    assert seq == "ATCGGATTAC"


def test_write_agp_creates_agp_creates_all_expected_rows_and_columns(fasta_split, tmp_path, write_fasta):
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", None), ("Y", "CC", None)])
    out = tmp_path / "out"

    fasta_split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=True,
        delete_existing_files=False,
        unique_file_names=False,
        force_max_seq_length=True,
        max_seqs_per_file=None,
        max_seq_length_per_file=4,
        min_chunk_length=None,
        max_files_per_directory=None,
        max_dirs_per_directory=None,
    )

    agp = out / "in.agp"
    assert agp.exists()

    lines = read_agp_lines(agp)
    assert lines[0].startswith("# AGP-version 2.0")

    comp_lines = [ln for ln in lines[1:] if ln and not ln.startswith("#")]
    assert len(comp_lines) == 4

    expected_values = [
        ["X", "1", "4", "1", "W", "X_chunk_start_0", "1", "4", "+"],
        ["X", "5", "8", "2", "W", "X_chunk_start_4", "1", "4", "+"],
        ["X", "9", "10", "3", "W", "X_chunk_start_8", "1", "2", "+"],
        ["Y", "1", "2", "1", "W", "Y", "1", "2", "+"],
    ]

    returned_values = [ln.split("\t") for ln in comp_lines]

    def agp_key(cols: list[str]) -> tuple[str, int]:
        return cols[0], int(cols[3])

    expected_values_sorted = sorted(expected_values, key=agp_key)
    returned_values_sorted = sorted(returned_values, key=agp_key)

    # Validate every column in every row
    assert returned_values_sorted == expected_values_sorted


def test_main_rejects_min_chunk_without_max_seq_length(fasta_split, tmp_path, write_fasta):
    in_fa = write_fasta("in.fa", [("a", "AA", None)])
    out = tmp_path / "out"

    # --min-chunk-length requires --max-seq-length-per-file
    with pytest.raises(ValueError, match="--min-chunk-length requires --max-seq-length-per-file"):
        fasta_split.parse_args(
            [
                "--fasta-file",
                str(in_fa),
                "--out-dir",
                str(out),
                "--min-chunk-length",
                "5",
            ]
        )


def test_parse_args_minimal_required(fasta_split, write_fasta):
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    args = fasta_split.parse_args(
        [
            "--fasta-file",
            str(in_fa),
        ]
    )

    assert args.fasta_file == in_fa
    assert args.write_agp is False
    assert args.force_max_seq_length is False


def test_parse_args_boolean_flags(fasta_split, write_fasta):
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    args = fasta_split.parse_args(
        [
            "--fasta-file",
            str(in_fa),
            "--write-agp",
            "--force-max-seq-length",
            "--unique-file-names",
        ]
    )

    assert args.write_agp
    assert args.force_max_seq_length
    assert args.unique_file_names


def test_parse_args_numeric_arguments(fasta_split, write_fasta):
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    args = fasta_split.parse_args(
        [
            "--fasta-file",
            str(in_fa),
            "--max-seqs-per-file",
            "5",
            "--max-seq-length-per-file",
            "100",
        ]
    )

    assert args.max_seqs_per_file == 5
    assert args.max_seq_length_per_file == 100


def test_parse_args_numeric_min_value_enforced(fasta_split, write_fasta):
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    with pytest.raises(SystemExit):
        fasta_split.parse_args(
            [
                "--fasta-file",
                str(in_fa),
                "--max-seqs-per-file",
                "0",
            ]
        )
