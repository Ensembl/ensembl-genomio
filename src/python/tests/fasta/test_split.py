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
import logging
from pathlib import Path
import re

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.io.genomio.fasta import split
from .utils import force_open_failure_for_suffix, read_fasta


def list_output_fastas(out_dir: Path) -> list[Path]:
    return sorted(p for p in out_dir.rglob("*.fa") if p.is_file())


def make_record(record_id: str, description: str) -> SeqRecord:
    r = SeqRecord(Seq("ACGT"), id=record_id)
    r.description = description
    return r


def read_agp_lines(agp_path: Path) -> list[str]:
    return agp_path.read_text(encoding="utf-8").splitlines()


@pytest.mark.parametrize(
    "name,expected",
    [("in.fa", "in"), ("in.fa.gz", "in"), ("in", "in")],
)
def test_get_fasta_basename(tmp_path: Path, name: str, expected: str) -> None:
    p = tmp_path / name
    p.write_text(">x\nA\n", encoding="utf-8")
    assert split._get_fasta_basename(p) == expected


def test_desc_without_id_header_only() -> None:
    rec = make_record("seq1", "seq1")
    assert split._description_without_id(rec) == ""


def test_desc_without_id_with_description() -> None:
    rec = make_record("seq1", "seq1 some annotation")
    assert split._description_without_id(rec) == "some annotation"


def test_desc_without_id_repeated_id_in_description() -> None:
    # original header: >seq1 seq1 description
    rec = make_record("seq1", "seq1 seq1 description")
    assert split._description_without_id(rec) == "seq1 description"


def test_desc_without_id_returns_original_when_id_not_prefix() -> None:
    rec = make_record("seq1", "other description")
    assert split._description_without_id(rec) == "other description"


def test_clean_previous_output_deletes_numeric_top_level_dirs(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("seq1", "ACGT", None)])
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    # Create numeric dirs with valid outputs
    d1 = out_dir / "1"
    d2 = out_dir / "2"
    d1.mkdir()
    d2.mkdir()
    fasta1 = d1 / "in.1.fa"
    fasta1.write_text(">x\nA\n", encoding="utf-8")
    fasta2 = d2 / "in.2.1.fa"
    fasta2.write_text(">y\nC\n", encoding="utf-8")

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

    assert fasta1.exists()
    assert fasta2.exists()

    split._clean_previous_output(in_fa, out_dir)

    assert not d1.exists()
    assert not d2.exists()
    assert other.exists()
    assert other_numeric.exists()
    assert other_zero.exists()
    assert (other / "in.3.fa").exists()
    assert (other_numeric / "in.4.fa").exists()
    assert (other_zero / "in.5.fa").exists()


def test_clean_previous_output_unlinks_agp_when_present(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("seq1", "ACGT", None)])
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    d1 = out_dir / "1"
    d1.mkdir()
    (d1 / "in.1.fa").write_text(">x\nA\n", encoding="utf-8")

    agp = out_dir / "in.agp"
    agp.write_text("# AGP-version 2.0\n", encoding="utf-8")
    assert agp.exists()

    split._clean_previous_output(in_fa, out_dir)

    assert not agp.exists()


def test_clean_previous_output_aborts_on_unexpected_file_and_deletes_nothing(
    tmp_path: Path, write_fasta
) -> None:
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
        split._clean_previous_output(in_fa, out_dir)

    # Nothing deleted because validation failed before deletion loop
    assert d1.exists()
    assert d2.exists()
    assert (d1 / "in.1.fa").exists()
    assert (d2 / "unexpected.fa").exists()


def test_clean_previous_output_returns_when_out_dir_missing(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out_dir = tmp_path / "does_not_exist"
    assert not out_dir.exists()

    split._clean_previous_output(in_fa, out_dir)


def test_check_contents_deletable_rejects_unexpected_nested_file(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    d1 = out_dir / "1"
    d2 = d1 / "1"
    d2.mkdir(parents=True)

    (d2 / "unexpected.txt").write_text("x", encoding="utf-8")

    output_file_re = re.compile(r"^in\.[1-9]\d*(\.[1-9]\d*)?\.fa$")
    with pytest.raises(RuntimeError, match=r"Unexpected file identified"):
        split._check_contents_deletable(d1, output_file_re)

    assert d2.exists()
    assert (d2 / "unexpected.txt").exists()


def test_outputwriter_basename_and_first_file_created(tmp_path: Path) -> None:
    in_fa = tmp_path / "in.fa.gz"
    in_fa.write_text(">x\nACGT\n", encoding="utf-8")

    out = tmp_path / "out"
    w = split.OutputWriter(
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
    write_fasta,
    tmp_path: Path,
    max_dirs_per_directory: int | None,
    dir_index: int,
    expected: list[str],
) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = split.OutputWriter(
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


def test_file_and_dir_index_with_max_files(write_fasta, tmp_path: Path) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = split.OutputWriter(
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


def test_unique_file_names_include_dir_index(write_fasta, tmp_path: Path) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = split.OutputWriter(
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


def test_add_agp_entry_when_agp_disabled(write_fasta, tmp_path: Path) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    writer = split.OutputWriter(
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


def test_write_record_requires_agp_args_when_write_agp_true(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("x", "ACGT", None)])
    out = tmp_path / "out"
    w = split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=True,
        unique_file_names=False,
    )
    try:
        rec = SeqRecord(Seq("ACGT"), id="x")
        with pytest.raises(ValueError, match="All AGP arguments must be provided"):
            w.write_record(rec, agp_object_id="x")  # missing agp_start/agp_end/agp_part_nr
    finally:
        w.close()


def test_split_fasta_empty_input_no_outputs(tmp_path: Path) -> None:
    in_fa = tmp_path / "empty.fa"
    in_fa.write_bytes(b"")
    out = tmp_path / "out"

    split.split_fasta(
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


def test_split_by_max_seqs_per_file_rollover(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta(
        "in.fa",
        [("a", "AA", None), ("b", "CC", None), ("c", "GG", None), ("d", "TT", None), ("e", "AT", None)],
    )
    out = tmp_path / "out"

    split.split_fasta(
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


def test_split_by_max_seq_length_per_file_rollover(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("a", "AAAA", None), ("b", "TTTT", None)])
    out = tmp_path / "out"

    split.split_fasta(
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


def test_force_chunking_splits_long_record(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", "desc")])
    out = tmp_path / "out"

    split.split_fasta(
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


def test_force_chunking_merges_small_remainder(
    tmp_path: Path, write_fasta, caplog: pytest.LogCaptureFixture
) -> None:
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", None)])
    out = tmp_path / "out"

    split.split_fasta(
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


def test_force_chunking_does_not_merge_when_remainder_meets_min_chunk_length(
    tmp_path: Path, write_fasta, caplog: pytest.LogCaptureFixture
) -> None:
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", None)])  # len=10
    out = tmp_path / "out"

    split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        force_max_seq_length=True,
        max_seq_length_per_file=4,
        min_chunk_length=2,  # last chunk length is 2 -> should NOT merge
    )

    assert not any("merging with previous chunk" in r.message for r in caplog.records)

    fastas = list_output_fastas(out)
    ids = []
    seqs = []
    for fasta in fastas:
        for record_id, seq in read_fasta(fasta).items():
            ids.append(record_id)
            seqs.append(seq)

    assert ids == ["X_chunk_start_0", "X_chunk_start_4", "X_chunk_start_8"]
    assert seqs == ["ATCG", "GATT", "AC"]


def test_overlong_record_without_chunking_warns_and_written_whole(
    tmp_path: Path, write_fasta, caplog: pytest.LogCaptureFixture
) -> None:
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", None)])
    out = tmp_path / "out"

    split.split_fasta(
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


def test_overlong_record_without_chunking_rolls_over_if_file_already_has_record(
    tmp_path: Path, write_fasta, caplog: pytest.LogCaptureFixture
) -> None:
    # first record short, second record overlong
    in_fa = write_fasta("in.fa", [("A", "AA", None), ("X", "ATCGGATTAC", None)])
    out = tmp_path / "out"

    split.split_fasta(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        force_max_seq_length=False,
        max_seq_length_per_file=4,
    )

    assert any("chunking not enabled" in r.message for r in caplog.records)

    fastas = list_output_fastas(out)
    assert len(fastas) == 2


def test_split_fasta_calls_clean_previous_output_when_requested(
    tmp_path: Path, write_fasta, monkeypatch: pytest.MonkeyPatch
) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    called = {"clean": False}

    def fake_clean(fasta_file: Path, out_dir_arg: Path) -> None:
        called["clean"] = True

    monkeypatch.setattr(split, "_clean_previous_output", fake_clean)

    split.split_fasta(
        fasta_file=in_fa,
        out_dir=out_dir,
        delete_existing_files=True,
    )

    assert called["clean"] is True


def test_create_agp_file_open_failure_raises_runtimeerror(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, write_fasta
) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    monkeypatch.setattr("builtins.open", force_open_failure_for_suffix(".agp"))

    with pytest.raises(RuntimeError, match=r"Failed to open AGP file"):
        split.OutputWriter(
            fasta_file=in_fa,
            out_dir=out,
            write_agp=True,
            unique_file_names=False,
        )


def test_create_agp_file_returns_if_agp_file_none(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    w = split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=False,
        unique_file_names=False,
    )
    try:
        w.create_agp_file()
        assert not (out / "in.agp").exists()
    finally:
        w.close()


def test_open_new_file_open_failure_raises_runtimeerror(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, write_fasta
) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    monkeypatch.setattr("builtins.open", force_open_failure_for_suffix(".fa"))

    with pytest.raises(RuntimeError, match=r"Failed to open output file"):
        split.OutputWriter(
            fasta_file=in_fa,
            out_dir=out,
            write_agp=False,
            unique_file_names=False,
        )


def test_write_agp_creates_agp_creates_all_expected_rows_and_columns(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("X", "ATCGGATTAC", None), ("Y", "CC", None)])
    out = tmp_path / "out"

    split.split_fasta(
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


def test_outputwriter_close_closes_open_handles(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("x", "A", None)])
    out = tmp_path / "out"

    w = split.OutputWriter(
        fasta_file=in_fa,
        out_dir=out,
        write_agp=True,
        unique_file_names=False,
    )

    assert w._fh is not None
    assert not w._fh.closed
    assert w._agp_fh is not None
    assert not w._agp_fh.closed

    w.close()
    assert w._fh is None
    assert w._agp_fh is None

    # Closing again should be no-op and not raise an exception
    w.close()
    assert w._fh is None
    assert w._agp_fh is None


def test_parse_args_rejects_min_chunk_without_max_seq_length(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("a", "AA", None)])
    out = tmp_path / "out"

    # --min-chunk-length requires --max-seq-length-per-file
    with pytest.raises(ValueError, match="--min-chunk-length requires --max-seq-length-per-file"):
        split.parse_args(
            [
                "--fasta-file",
                str(in_fa),
                "--out-dir",
                str(out),
                "--min-chunk-length",
                "5",
            ]
        )


def test_parse_args_minimal_required(write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    args = split.parse_args(
        [
            "--fasta-file",
            str(in_fa),
        ]
    )

    assert args.fasta_file == in_fa
    assert args.write_agp is False
    assert args.force_max_seq_length is False


def test_parse_args_boolean_flags(write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    args = split.parse_args(
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


def test_parse_args_numeric_arguments(write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    args = split.parse_args(
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


def test_parse_args_numeric_min_value_enforced(write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("a", "AA", None)])

    with pytest.raises(SystemExit):
        split.parse_args(
            [
                "--fasta-file",
                str(in_fa),
                "--max-seqs-per-file",
                "0",
            ]
        )


def test_parse_args_out_dir_is_path(tmp_path: Path, write_fasta) -> None:
    in_fa = write_fasta("in.fa", [("a", "AA", None)])
    out = tmp_path / "out"

    args = split.parse_args(["--fasta-file", str(in_fa), "--out-dir", str(out)])
    assert isinstance(args.out_dir, Path)
    assert args.out_dir == out


def test_main_logs_and_reraises_exceptions(
    monkeypatch: pytest.MonkeyPatch, write_fasta, caplog: pytest.LogCaptureFixture
) -> None:
    in_fa = write_fasta("in.fa", [("a", "AA", None)])
    monkeypatch.setattr(split, "init_logging_with_args", lambda args: None)

    def raise_main_exception(*args: object, **kwargs: object) -> None:
        raise RuntimeError("Simulated exception in main")

    monkeypatch.setattr(split, "split_fasta", raise_main_exception)

    caplog.set_level(logging.ERROR, logger="root")

    with pytest.raises(RuntimeError, match="Simulated exception in main"):
        split.main(["--fasta-file", str(in_fa)])
    assert "Error processing FASTA file" in caplog.text
