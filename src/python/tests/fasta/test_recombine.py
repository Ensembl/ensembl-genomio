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

import importlib
from pathlib import Path
import pytest
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.io.genomio.fasta.recombine import main as FastaRecombine


@pytest.fixture(scope="session")
def fasta_recombine():
    return importlib.import_module("ensembl.io.genomio.fasta.recombine")


@pytest.fixture
def chunk_re():
    return re.compile(r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$")


def read_fasta(path: Path) -> dict[str, str]:
    with open(path, "r", encoding="utf-8") as fh:
        return {r.id: str(r.seq) for r in SeqIO.parse(fh, "fasta")}


def write_manifest(path: Path, lines: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_get_fasta_paths_preserves_manifest_order(fasta_recombine, write_fasta, tmp_path):
    a = write_fasta("a.fa", [("r1", "AAA", None)])
    b = write_fasta("b.fa", [("r2", "BBB", None)])

    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            str(b),
            str(a),
        ],
    )

    paths = fasta_recombine._get_fasta_paths(manifest)
    assert paths == [b.resolve(), a.resolve()]


def test_get_fasta_paths_ignores_comments_and_blank_lines(fasta_recombine, write_fasta, tmp_path):
    a = write_fasta("a.fa", [("r1", "AAA", None)])
    b = write_fasta("b.fa.gz", [("r2", "TT", None)], gz=True)

    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            "# comment",
            "",
            str(a),
            "   ",
            str(b),
            "# trailing comment",
        ],
    )

    paths = fasta_recombine._get_fasta_paths(manifest)
    assert paths == [a.resolve(), b.resolve()]


def test_get_fasta_paths_relative_paths_are_anchored_to_manifest_dir(fasta_recombine, write_fasta, tmp_path):
    sub = tmp_path / "sub"
    a = write_fasta("sub/a.fa", [("r1", "AAA", None)])

    manifest = write_manifest(
        sub / "manifest.txt",
        [
            "a.fa",  # relative to manifest dir
        ],
    )

    paths = fasta_recombine._get_fasta_paths(manifest)
    assert paths == [a.resolve()]


def test_get_fasta_paths_missing_entry_raises_with_line_number(fasta_recombine, write_fasta, tmp_path):
    a = write_fasta("a.fa", [("r1", "AAA", None)])
    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            str(a),
            "missing.fa",
        ],
    )

    with pytest.raises(FileNotFoundError, match=r"line 2"):
        fasta_recombine._get_fasta_paths(manifest)


def test_get_fasta_paths_directory_entry_raises(fasta_recombine, tmp_path):
    d = tmp_path / "adir"
    d.mkdir()

    manifest = write_manifest(
        tmp_path / "manifest.txt",
        [
            str(d),
        ],
    )

    with pytest.raises(ValueError, match=r"Manifest entry.*line 1"):
        fasta_recombine._get_fasta_paths(manifest)


def test_get_fasta_paths_manifest_must_be_a_file(fasta_recombine, tmp_path):
    manifest_dir = tmp_path / "manifest_dir"
    manifest_dir.mkdir()

    with pytest.raises(ValueError, match="Manifest is not a file"):
        fasta_recombine._get_fasta_paths(manifest_dir)


def test_fasta_record_cache_load_and_switch_files(fasta_recombine, write_fasta):
    f1 = write_fasta("one.fa", [("a", "AAAA", None), ("b", "CC", None)])
    f2 = write_fasta("two.fa", [("c", "GGG", None)])

    loc1a = fasta_recombine.RecordLocation(path=f1, description="a desc")
    loc2c = fasta_recombine.RecordLocation(path=f2, description="c desc")

    cache = fasta_recombine.FastaRecordCache()
    rec_a = cache.get("a", loc1a)
    assert rec_a.id == "a"
    assert str(rec_a.seq) == "AAAA"

    # switch file
    rec_c = cache.get("c", loc2c)
    assert rec_c.id == "c"
    assert str(rec_c.seq) == "GGG"


def test_fasta_record_cache_missing_record_raises_keyerror(fasta_recombine, write_fasta):
    f = write_fasta("one.fa", [("a", "AAAA", None)])
    loc = fasta_recombine.RecordLocation(path=f, description="x")
    cache = fasta_recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="Record 'invalid' not found"):
        cache.get("invalid", loc)


def test_fasta_record_cache_duplicate_id_within_file_raises(fasta_recombine, write_fasta):
    f = write_fasta("dup.fa", [("a", "AAAA", None), ("a", "TTTT", None)])
    cache = fasta_recombine.FastaRecordCache()
    with pytest.raises(ValueError, match="Duplicate record id 'a'"):
        cache._load_file(f)


def test_build_index_collects_locations_first_seen_and_chunks(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta(
        "d/f1.fa",
        [
            ("X_chunk_start_0", "AAAA", "chunk0"),
            ("Y", "CC", "unchunked"),
        ],
    )
    f2 = write_fasta(
        "d/f2.fa",
        [
            ("X_chunk_start_4", "TT", "chunk4"),
            ("Z", "GG", "z"),
        ],
    )

    locations, first_seen, chunks = fasta_recombine._build_index(chunk_re, [f1, f2])

    assert "Y" in locations
    assert "X_chunk_start_0" in locations
    assert first_seen["X"] < first_seen["Y"]  # X base first observed in f1
    assert chunks["X"] == [(0, "X_chunk_start_0"), (4, "X_chunk_start_4")]


def test_build_index_duplicate_record_id_raises(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta("f1.fa", [("A", "AAA", None)])
    f2 = write_fasta("f2.fa", [("A", "TTT", None)])  # duplicate id across files
    with pytest.raises(ValueError, match="Duplicate FASTA record id encountered during indexing"):
        fasta_recombine._build_index(chunk_re, [f1, f2])


def test_build_index_no_headers_raises(fasta_recombine, chunk_re, tmp_path):
    # create an empty file that parses to no records
    p = tmp_path / "empty.fa"
    p.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="No FASTA headers found"):
        fasta_recombine._build_index(chunk_re, [p])


def test_records_from_headers_passes_unchunked(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta("f1.fa", [("A", "AAA", "descA")])
    locations, first_seen, chunks = fasta_recombine._build_index(chunk_re, [f1])

    cache = fasta_recombine.FastaRecordCache()
    recs = list(fasta_recombine._records_from_headers(locations, first_seen, chunks, cache))
    assert [(r.id, str(r.seq)) for r in recs] == [("A", "AAA")]
    assert recs[0].description.endswith("descA")


def test_records_from_headers_reassembles_contiguous_chunks(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta(
        "f1.fa",
        [
            ("X_chunk_start_0", "AAAA", "c0"),
            ("X_chunk_start_4", "TT", "c4"),
        ],
    )
    locations, first_seen, chunks = fasta_recombine._build_index(chunk_re, [f1])
    cache = fasta_recombine.FastaRecordCache()

    recs = list(fasta_recombine._records_from_headers(locations, first_seen, chunks, cache))
    assert len(recs) == 1
    assert recs[0].id == "X"
    assert str(recs[0].seq) == "AAAATT"


def test_records_from_headers_noncontiguous_chunks_raises(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta(
        "f1.fa",
        [
            ("X_chunk_start_0", "AAAA", None),
            ("X_chunk_start_99", "TT", None),  # wrong start
        ],
    )
    locations, first_seen, chunks = fasta_recombine._build_index(chunk_re, [f1])
    cache = fasta_recombine.FastaRecordCache()

    with pytest.raises(ValueError, match="Non-contiguous chunks for 'X'"):
        list(fasta_recombine._records_from_headers(locations, first_seen, chunks, cache))


def test_records_from_headers_missing_base_record_raises(fasta_recombine):
    locations = {}
    first_seen = {"A": 0}
    chunks = {}
    cache = fasta_recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="Base record 'A' not found"):
        list(fasta_recombine._records_from_headers(locations, first_seen, chunks, cache))


def test_records_from_headers_missing_chunk_record_raises(fasta_recombine):
    # base "X" has chunks but one chunk id isn't in locations
    locations = {}
    first_seen = {"X": 0}
    chunks = {"X": [(0, "X_chunk_start_0")]}
    cache = fasta_recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="Chunk record 'X_chunk_start_0' not found"):
        list(fasta_recombine._records_from_headers(locations, first_seen, chunks, cache))


def test_parse_agp_ignores_comments_and_blank_lines(fasta_recombine, tmp_path):
    agp = tmp_path / "x.agp"
    agp.write_text(
        "# comment\n\n" "obj1\t1\t4\t1\tW\tpart1\t1\t4\t+\n",
        encoding="utf-8",
    )
    out = fasta_recombine._parse_agp(agp, allow_revcomp=False)
    assert "obj1" in out
    assert out["obj1"][0].part_id == "part1"


def test_parse_agp_rejects_short_line(fasta_recombine, tmp_path):
    agp = tmp_path / "bad.agp"
    agp.write_text("obj\t1\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected >= 9 columns"):
        fasta_recombine._parse_agp(agp, allow_revcomp=False)


def test_parse_agp_rejects_non_W_component(fasta_recombine, tmp_path):
    agp = tmp_path / "bad.agp"
    agp.write_text("obj\t1\t2\t1\tN\tgap\t1\t2\t+\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported AGP component type"):
        fasta_recombine._parse_agp(agp, allow_revcomp=False)


def test_parse_agp_rejects_minus_orientation_unless_allowed(fasta_recombine, tmp_path):
    agp = tmp_path / "bad.agp"
    agp.write_text("obj\t1\t2\t1\tW\tpart\t1\t2\t-\n", encoding="utf-8")
    with pytest.raises(ValueError, match="--allow-revcomp is not enabled"):
        fasta_recombine._parse_agp(agp, allow_revcomp=False)

    ok = fasta_recombine._parse_agp(agp, allow_revcomp=True)
    assert ok["obj"][0].orientation == "-"


def test_parse_agp_empty_file_raises(fasta_recombine, tmp_path):
    agp = tmp_path / "empty.agp"
    agp.write_text("# only comment\n", encoding="utf-8")
    with pytest.raises(ValueError, match="contained no component lines"):
        fasta_recombine._parse_agp(agp, allow_revcomp=False)


def test_agp_component_seq_plus_and_minus(fasta_recombine):
    record = SeqRecord(Seq("ACGTGGTT"), "part", "test_record")
    assert str(fasta_recombine._agp_component_seq(record, 1, 4, "+", allow_revcomp=False)) == "ACGT"
    assert str(fasta_recombine._agp_component_seq(record, 3, 6, "-", allow_revcomp=True)) == str(
        Seq("GTGG").reverse_complement()
    )

    with pytest.raises(ValueError, match="--allow-revcomp is not enabled"):
        fasta_recombine._agp_component_seq(record, 3, 6, "-", allow_revcomp=False)

    with pytest.raises(ValueError, match="Invalid AGP orientation"):
        fasta_recombine._agp_component_seq(record, 1, 2, "?", allow_revcomp=True)


def test_records_from_agp_success(fasta_recombine, write_fasta, chunk_re, tmp_path):
    parts = write_fasta("parts.fa", [("p1", "AAAA", "p1desc"), ("p2", "TT", "p2desc")])
    locations, _, _ = fasta_recombine._build_index(chunk_re, [parts])

    # AGP describes obj = p1(1..4) + p2(1..2)
    agp = tmp_path / "x.agp"
    agp.write_text(
        "obj\t1\t4\t1\tW\tp1\t1\t4\t+\n" "obj\t5\t6\t2\tW\tp2\t1\t2\t+\n",
        encoding="utf-8",
    )
    agp_entries = fasta_recombine._parse_agp(agp, allow_revcomp=False)

    cache = fasta_recombine.FastaRecordCache()
    recs = list(fasta_recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))
    assert [(r.id, str(r.seq)) for r in recs] == [("obj", "AAAATT")]


def test_records_from_agp_sorts_parts_out_of_order(fasta_recombine, write_fasta, chunk_re, tmp_path):
    parts = write_fasta("parts.fa", [("p1", "AAAA", "p1desc"), ("p2", "TT", "p2desc")])
    locations, _, _ = fasta_recombine._build_index(chunk_re, [parts])

    # AGP lines deliberately out-of-order (p2 line first); should still assemble p1 + p2
    agp = tmp_path / "x.agp"
    agp.write_text(
        "obj\t5\t6\t2\tW\tp2\t1\t2\t+\n" "obj\t1\t4\t1\tW\tp1\t1\t4\t+\n",
        encoding="utf-8",
    )
    agp_entries = fasta_recombine._parse_agp(agp, allow_revcomp=False)

    cache = fasta_recombine.FastaRecordCache()
    recs = list(fasta_recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))
    assert [(r.id, str(r.seq)) for r in recs] == [("obj", "AAAATT")]


def test_records_from_agp_missing_component_raises(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta("parts.fa", [("p1", "AAAA", None)])
    locations, _, _ = fasta_recombine._build_index(chunk_re, [f1])

    agp_entries = {
        "obj": [
            fasta_recombine.AgpEntry(
                record="obj",
                record_start=1,
                record_end=4,
                part_number=1,
                part_id="MISSING",
                part_start=1,
                part_end=4,
                orientation="+",
            )
        ]
    }

    cache = fasta_recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="not found in indexed FASTA headers"):
        list(fasta_recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))


def test_records_from_agp_noncontiguous_object_raises(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta("parts.fa", [("p1", "AAAA", None), ("p2", "TT", None)])
    locations, _, _ = fasta_recombine._build_index(chunk_re, [f1])

    agp_entries = {
        "obj": [
            fasta_recombine.AgpEntry("obj", 1, 4, 1, "p1", 1, 4, "+"),
            # should start at 5, but starts at 6
            fasta_recombine.AgpEntry("obj", 6, 7, 2, "p2", 1, 2, "+"),
        ]
    }
    cache = fasta_recombine.FastaRecordCache()
    with pytest.raises(ValueError, match="Non-contiguous AGP"):
        list(fasta_recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))


def test_records_from_agp_length_mismatch_raises(fasta_recombine, write_fasta, chunk_re):
    f1 = write_fasta("parts.fa", [("p1", "AAAA", None)])
    locations, _, _ = fasta_recombine._build_index(chunk_re, [f1])

    # asks for 1..10 but record shorter => extracted length != expected_len
    agp_entries = {"obj": [fasta_recombine.AgpEntry("obj", 1, 10, 1, "p1", 1, 10, "+")]}
    cache = fasta_recombine.FastaRecordCache()
    with pytest.raises(ValueError, match="Length mismatch"):
        list(fasta_recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))


def test_recombine_fasta_header_mode_end_to_end_manifest_order(
    fasta_recombine, write_fasta, chunk_re, tmp_path
):
    f1 = write_fasta("f1.fa", [("Y", "CC", None)])
    f2 = write_fasta("f2.fa", [("X_chunk_start_0", "AAAA", None), ("X_chunk_start_4", "TT", None)])

    manifest = write_manifest(tmp_path / "manifest.txt", [str(f1), str(f2)])

    out = tmp_path / "out.fa"
    fasta_recombine.recombine_fasta(
        fasta_manifest=manifest,
        out_fasta=out,
        chunk_re=chunk_re,
        agp_file=None,
        allow_revcomp=False,
    )

    seqs = read_fasta(out)
    assert seqs["X"] == "AAAATT"
    assert seqs["Y"] == "CC"

    ids_in_file = [r.id for r in SeqIO.parse(out, "fasta")]
    assert ids_in_file == ["Y", "X"]


def test_recombine_fasta_agp_mode_end_to_end(fasta_recombine, write_fasta, chunk_re, tmp_path):
    parts = write_fasta("parts.fa", [("p1", "AACCGG", None), ("p2", "TTAA", None)])

    agp = tmp_path / "x.agp"
    agp.write_text(
        "obj\t1\t2\t1\tW\tp1\t1\t2\t+\n" "obj\t3\t4\t2\tW\tp2\t3\t4\t-\n",  # part2 slice 'AA' revcomp => 'TT'
        encoding="utf-8",
    )

    manifest = write_manifest(tmp_path / "manifest.txt", [str(parts)])

    out = tmp_path / "out.fa"
    fasta_recombine.recombine_fasta(
        fasta_manifest=manifest,
        out_fasta=out,
        chunk_re=chunk_re,
        agp_file=agp,
        allow_revcomp=True,
    )
    seqs = read_fasta(out)
    assert seqs["obj"] == "AATT"


def test_validate_regex_accepts_alternative_pattern(fasta_recombine):
    # Alternative chunk naming scheme: "<base>.part.<start>_<end>"
    chunk_re = fasta_recombine._validate_regex(r"^(?P<base>.+)\.part\.(?P<start>\d+)_(?P<end>\d+)$")
    m = chunk_re.match("contigA.part.11_20")
    assert m is not None
    assert m.group("base") == "contigA"
    assert m.group("start") == "11"


def test_validate_regex_rejects_invalid_pattern(fasta_recombine):
    with pytest.raises(ValueError, match="Invalid --chunk-id-regex"):
        fasta_recombine._validate_regex(r"^(?P<base>.+)_(?P<start>\d+$")  # missing ')'


def test_validate_regex_requires_named_groups(fasta_recombine):
    with pytest.raises(ValueError, match="must define named capture groups"):
        fasta_recombine._validate_regex(r"^(.+)_(\d+)$")

    ok = fasta_recombine._validate_regex(r"^(?P<base>.+)_(?P<start>\d+)$")
    assert ok.match("X_12")


def test_recombine_header_mode_with_alternative_regex_end_to_end(
    fasta_recombine, write_fasta, tmp_path: Path
):
    chunk_re = fasta_recombine._validate_regex(r"^(?P<base>.+)\.chunk\.(?P<start>\d+)$")

    f1 = write_fasta("d/a.fa", [("X.chunk.0", "AAAA", None)])
    f2 = write_fasta("d/b.fa", [("X.chunk.4", "TT", None)])
    f3 = write_fasta("d/c.fa", [("Y", "CC", None)])

    manifest = write_manifest(tmp_path / "manifest.txt", [str(f1), str(f2), str(f3)])

    out = tmp_path / "out.fa"
    fasta_recombine.recombine_fasta(
        fasta_manifest=manifest,
        out_fasta=out,
        chunk_re=chunk_re,
        agp_file=None,
        allow_revcomp=False,
    )

    seqs = read_fasta(out)
    assert seqs["X"] == "AAAATT"
    assert seqs["Y"] == "CC"


def test_parse_args_minimal_required(fasta_recombine, write_fasta, tmp_path):
    f = write_fasta("a.fa", [("A", "AAA", None)])
    manifest = write_manifest(tmp_path / "manifest.txt", [str(f)])
    out = tmp_path / "out.fa"
    args = fasta_recombine.parse_args(
        [
            "--fasta-manifest",
            str(manifest),
            "--out-fasta",
            str(out),
        ]
    )

    assert args.fasta_manifest == manifest
    assert args.out_fasta == out

    assert args.allow_revcomp is False
    assert isinstance(args.chunk_id_regex, str)
    assert fasta_recombine._validate_regex(args.chunk_id_regex).match("X_chunk_start_0")


def test_parse_args_boolean_flags_manifest(fasta_recombine, write_fasta, tmp_path):
    f = write_fasta("a.fa", [("A", "AAA", None)])
    manifest = write_manifest(tmp_path / "manifest.txt", [str(f)])
    out = tmp_path / "out.fa"

    args = fasta_recombine.parse_args(
        [
            "--fasta-manifest",
            str(manifest),
            "--out-fasta",
            str(out),
            "--allow-revcomp",
        ]
    )
    assert args.allow_revcomp is True
