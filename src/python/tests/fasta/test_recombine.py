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
"""Unit testing of `ensembl.io.genomio.fasta.recombine` module."""

from contextlib import nullcontext as does_not_raise
from pathlib import Path
import pytest
from pytest import param
import re
from typing import ContextManager
from unittest.mock import Mock, patch


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from deepdiff import DeepDiff

from ensembl.io.genomio.fasta import recombine
from ensembl.io.genomio.utils.agp_utils import parse_agp


CHUNK_RE = re.compile(r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$")


def test_fasta_record_cache_load_and_switch_files(data_dir: Path) -> None:
    """
    Tests loading and switching of files by `recombine.FastaRecordCache()`.

    Args:
        data_dir: Module's test data directory fixture.
    """
    fasta_file1 = data_dir / "input_a.fa"
    fasta_file2 = data_dir / "input_b.fa"

    loc1 = recombine.RecordLocation(path=fasta_file1, description="a desc")
    loc2 = recombine.RecordLocation(path=fasta_file2, description="b desc")

    cache = recombine.FastaRecordCache()
    rec_a_seq1 = cache.get("a_seq1", loc1)
    assert rec_a_seq1.id == "a_seq1"
    assert str(rec_a_seq1.seq) == "AAAA"

    # switch file
    rec_b_seq1 = cache.get("b_seq1", loc2)
    assert rec_b_seq1.id == "b_seq1"
    assert str(rec_b_seq1.seq) == "GGG"


def test_fasta_record_cache_missing_record_raises_keyerror(data_dir: Path) -> None:
    """
    Tests `recombine.FastaRecordCache()` raises error when invalid record requested.

    Args:
        data_dir: Module's test data directory fixture.
    """
    fasta_file = data_dir / "input_a.fa"
    loc = recombine.RecordLocation(path=fasta_file, description="description")
    cache = recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="Record 'missing' not found"):
        cache.get("missing", loc)


def test_fasta_record_cache_duplicate_id_within_file_raises(data_dir: Path) -> None:
    """
    Tests `recombine.FastaRecordCache()` raises error when loading file with duplicate record ID.

    Args:
        data_dir: Module's test data directory fixture.
    """
    fasta_file = data_dir / "duplicate_id.fa"
    cache = recombine.FastaRecordCache()
    with pytest.raises(ValueError, match="Duplicate record id 'seq_a'"):
        cache._load_file(fasta_file)


def test_build_index_collects_locations_first_seen_and_chunks(data_dir: Path) -> None:
    """
    Tests `recombine._build_index()`, ensuring sequence order captured correctly.

    Args:
        data_dir: Module's test data directory fixture.
    """
    fasta_file1 = data_dir / "build_index" / "input_a.fa"
    fasta_file2 = data_dir / "build_index" / "input_b.fa"
    locations, first_seen, chunks = recombine._build_index(CHUNK_RE, [fasta_file1, fasta_file2])

    assert {"X_chunk_start_0", "X_chunk_start_4", "Y", "Z"} <= set(locations)
    assert first_seen["X"] < first_seen["Y"] < first_seen["Z"]
    assert chunks["X"] == [(0, "X_chunk_start_0"), (4, "X_chunk_start_4")]


def test_build_index_duplicate_record_id_raises(data_dir: Path) -> None:
    """
    Tests `recombine._build_index()` raises error when duplicate record IDs found in files being indexed.

    Args:
        data_dir: Module's test data directory fixture.
    """
    fasta_file1 = data_dir / "duplicate_across_files" / "duplicate_a.fa"
    fasta_file2 = data_dir / "duplicate_across_files" / "duplicate_b.fa"
    with pytest.raises(ValueError, match="Duplicate FASTA record id encountered during indexing"):
        recombine._build_index(CHUNK_RE, [fasta_file1, fasta_file2])


def test_build_index_no_headers_raises(tmp_path: Path) -> None:
    """
    Tests `recombine._build_index()` raises error when attempting to index empty files.

    Args:
        tmp_path: Temporary directory provided by pytest.
    """
    empty_file = tmp_path / "empty.fa"
    empty_file.touch()
    with pytest.raises(ValueError, match="No FASTA headers found"):
        recombine._build_index(CHUNK_RE, [empty_file])


@pytest.mark.parametrize(
    "input_filename,regex,expectation",
    [
        param(
            "0_based_start.fa",
            CHUNK_RE,
            does_not_raise(("X", "AAAATT")),
            id="reassembles_0_based_contiguous_chunks",
        ),
        param(
            "1_based_start.fa",
            CHUNK_RE,
            does_not_raise(("Y", "CCCCGG")),
            id="reassembles_1_based_contiguous_chunks",
        ),
        param(
            "alternative_regex.fa",
            re.compile(r"^(?P<base>.+)\.chunk\.(?P<start>\d+)$"),
            does_not_raise(("Z", "TCCTAGA")),
            id="reassembles_with_alternative_regex",
        ),
        param(
            "unchunked.fa",
            CHUNK_RE,
            does_not_raise(("A", "CATTGA")),
            id="unchunked_passes_through",
        ),
        param(
            "non_contiguous.fa",
            CHUNK_RE,
            pytest.raises(ValueError, match="Non-contiguous chunks for 'X'"),
            id="rejects_non_contiguous_chunks",
        ),
    ],
)
def test_recombine_records_from_headers(
    data_dir: Path,
    input_filename: str,
    regex: re.Pattern[str],
    expectation: ContextManager,
) -> None:
    """
    Tests the `recombine._records_from_headers()` function.

    Args:
        data_dir: Module's test data directory fixture.
        input_filename: Name of input FASTA file.
        regex: re.Pattern for matching base and start groups in FASTA headers.
        expectation: Context manager for the expected outcome of the test (exception or not).
    """
    input_fa = data_dir / "from_headers" / input_filename
    with expectation as expected:
        locations, first_seen, chunks = recombine._build_index(regex, [input_fa])
        cache = recombine.FastaRecordCache()

        records = list(recombine._records_from_headers(locations, first_seen, chunks, cache))
        record_id, record_seq = expected
        assert len(records) == 1
        assert records[0].id == record_id
        assert str(records[0].seq) == record_seq


@pytest.mark.parametrize(
    "first_seen,chunks,key_error_str",
    [
        param(
            {"A": 0},
            {},
            "Base record 'A' not found",
            id="missing_base_record_raises",
        ),
        param(
            {"X": 0},
            {"X": [(0, "X_chunk_start_0")]},
            "Chunk record 'X_chunk_start_0' not found",
            id="missing_chunk_record_raises",
        ),
    ],
)
def test_records_from_headers_exceptions(
    first_seen: dict[str, int],
    chunks: dict[str, list[tuple[int, str]]],
    key_error_str: str,
) -> None:
    """
    Tests exceptions raised by the recombine._records_from_headers() function.

    Args:
        first_seen: dictionary representing order in which records observed.
        chunks: mapping of base record ID to a list of (start, chunk_record_id) pairs.
        key_error: expected string matching raised KeyError.
    """
    locations = {}
    cache = recombine.FastaRecordCache()
    with pytest.raises(KeyError, match=key_error_str):
        list(recombine._records_from_headers(locations, first_seen, chunks, cache))


@pytest.mark.parametrize(
    "orientation,allow_revcomp,expectation",
    [
        param(
            "+",
            False,
            does_not_raise("GTGG"),
            id="retrieves_forward_sequence",
        ),
        param("-", True, does_not_raise("CCAC"), id="retrieves_reverse_complement_sequence"),
        param(
            "-",
            False,
            pytest.raises(ValueError, match="--allow-revcomp is not enabled"),
            id="rejects_minus_strand_when_not_allowed",
        ),
        param(
            "?",
            True,
            pytest.raises(ValueError, match="Invalid AGP orientation"),
            id="rejects_invalid_orientation",
        ),
    ],
)
def test_agp_component_seq(
    orientation: str,
    allow_revcomp: bool,
    expectation: ContextManager,
) -> None:
    """
    Tests the `recombine._agp_component_seq()` function.

    Args:
        orientation: '+' or '-' orientation from AGP.
        allow_revcomp: Boolean indicating whether minus strand entries are permitted.
        expectation: Context manager for the expected outcome of the test (exception or not).
    """
    record = SeqRecord(Seq("ACGTGGTT"), "part", "test_record")
    with expectation as expected:
        sequence: Seq = recombine._agp_component_seq(record, 3, 6, orientation, allow_revcomp)
        assert str(sequence) == expected


@pytest.mark.parametrize(
    "agp_filename,expectation",
    [
        param(
            "ordered.agp",
            does_not_raise(("obj", "AAAATT")),
            id="reassembles_from_ordered_agp",
        ),
        param(
            "unordered.agp",
            does_not_raise(("obj", "AAAATT")),
            id="reassembles_from_unordered_agp",
        ),
        param(
            "length_mismatch.agp",
            pytest.raises(ValueError, match="Length mismatch"),
            id="rejects_agp_with_length_mismatch",
        ),
        param(
            "missing.agp",
            pytest.raises(KeyError, match="not found in indexed FASTA headers"),
            id="rejects_agp_with_unknown_sequence",
        ),
        param(
            "non_contiguous.agp",
            pytest.raises(ValueError, match="Non-contiguous AGP"),
            id="rejects_agp_with_coordinate_gap",
        ),
    ],
)
def test_records_from_agp(
    data_dir: Path,
    agp_filename: str,
    expectation: ContextManager,
) -> None:
    """
    Tests the `recombine._records_from_agp()` function.

    Args:
        data_dir: Module's test data directory fixture.
        agp_filename: Name of input AGP file.
        expectation: Context manager for the expected outcome of the test (exception or not).
    """
    locations, _, _ = recombine._build_index(CHUNK_RE, [data_dir / "from_agp" / "input.fa"])
    agp_entries = parse_agp(data_dir / "from_agp" / agp_filename, allow_revcomp=False)
    cache = recombine.FastaRecordCache()
    with expectation as expected:
        records = list(recombine._records_from_agp(agp_entries, locations, cache, False))
        record_id, record_seq = expected
        assert records[0].id == record_id
        assert str(records[0].seq) == record_seq


@pytest.mark.parametrize(
    "test_dir_name,agp_filename,expected",
    [
        param(
            "from_agp",
            "ordered.agp",
            [{"id": "obj", "seq": "AAAATT"}],
            id="agp_driven_recombination",
        ),
        param(
            "from_headers",
            None,
            [{"id": "Y", "seq": "CCCCGG"}],
            id="header_driven_recombination",
        ),
    ],
)
def test_recombine_fasta(
    data_dir: Path,
    tmp_path: Path,
    test_dir_name: str,
    agp_filename: str | None,
    expected: list[dict[str, str]],
) -> None:
    """
    Tests the `recombine.recombine_fasta()` function.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        test_dir_name: Name of data_dir subfolder containing test files.
        agp_filename: Name of (optional) AGP file.
        expected: results of reading FASTA output into list of dictionaries keyed by 'id' and 'seq'.
    """
    manifest = data_dir / test_dir_name / "manifest.txt"
    out = tmp_path / "out.fa"
    agp_file = data_dir / test_dir_name / agp_filename if agp_filename else None
    recombine.recombine_fasta(
        fasta_manifest=manifest,
        out_fasta=out,
        chunk_re=CHUNK_RE,
        agp_file=agp_file,
        allow_revcomp=False,
    )

    assert [{"id": r.id, "seq": str(r.seq)} for r in SeqIO.parse(out, "fasta")] == expected


@pytest.mark.parametrize(
    "arg_list, expected_params",
    [
        param(
            ["--fasta-manifest", __file__, "--out-fasta", "out.fa"],
            {
                "fasta_manifest": str(__file__),
                "out_fasta": "out.fa",
                "allow_revcomp": False,
                "chunk_id_regex": CHUNK_RE.pattern,
                "log_level": "WARNING",
            },
            id="Default args",
        ),
        param(
            [
                "--fasta-manifest",
                __file__,
                "--out-fasta",
                "out.fa",
                "--agp-file",
                __file__,
                "--allow-revcomp",
                "--chunk-id-regex",
                r"^(?P<base>.+)_start_(?P<start>\d+)$",
            ],
            {
                "fasta_manifest": str(__file__),
                "out_fasta": "out.fa",
                "agp_file": str(__file__),
                "allow_revcomp": True,
                "chunk_id_regex": r"^(?P<base>.+)_start_(?P<start>\d+)$",
                "log_level": "WARNING",
            },
            id="New arg values",
        ),
    ],
)
def test_parse_args(arg_list: list[str], expected_params: dict[str, str | int | bool]) -> None:
    """
    Tests the `recombine.parse_args()` function.

    Args:
        arg_list: List of command line arguments to parse.
        expected_params: Expected dictionary of parameters from parsing arguments.
    """
    args = recombine.parse_args(arg_list)
    # DeepDiff is not able to compare two objects of Path type - need to convert them to string
    setattr(args, "fasta_manifest", str(args.fasta_manifest))
    setattr(args, "out_fasta", str(args.out_fasta))
    if hasattr(args, "agp_file"):
        setattr(args, "agp_file", str(args.agp_file))
    assert not DeepDiff(vars(args), expected_params)


@patch("ensembl.io.genomio.fasta.recombine.recombine_fasta")
def test_main(mock_recombine_fasta: Mock, tmp_path: Path) -> None:
    """
    Tests the `recombine.main()` function (entry point).

    Args:
        mock_recombine_fasta: Mock object for the `recombine.recombine_fasta()` function.
        tmp_path: Temporary directory provided by pytest.
    """
    fasta_path = tmp_path / "in.fa"
    manifest_path = tmp_path / "manifest.txt"
    manifest_path.touch()
    recombine.main(["--fasta-manifest", str(manifest_path), "--out-fasta", str(fasta_path)])
    mock_recombine_fasta.assert_called_once_with(
        fasta_manifest=manifest_path,
        out_fasta=fasta_path,
        chunk_re=re.compile(r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$"),
        agp_file=None,
        allow_revcomp=False,
    )


@patch("ensembl.io.genomio.fasta.recombine.recombine_fasta")
def test_main_raises_exception(mock_recombine_fasta: Mock, tmp_path: Path) -> None:
    """
    Tests the `recombine.main()` function (entry point) raises exception on error.

    Args:
        mock_recombine_fasta: Mock object for the `recombine.recombine_fasta()` function.
        tmp_path: Temporary directory provided by pytest.
    """
    fasta_path = tmp_path / "in.fa"
    manifest_path = tmp_path / "manifest.txt"
    manifest_path.touch()
    mock_recombine_fasta.side_effect = RuntimeError("Mocked exception")
    with pytest.raises(RuntimeError, match=r"Mocked exception"):
        recombine.main(["--fasta-manifest", str(manifest_path), "--out-fasta", str(fasta_path)])
