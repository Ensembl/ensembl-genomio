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

from pathlib import Path
import pytest
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.io.genomio.fasta import recombine
from ensembl.io.genomio.utils.agp_utils import parse_agp
from ensembl.io.genomio.utils.chunk_utils import validate_regex

from .._helpers import write_manifest


CHUNK_RE = re.compile(r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$")


def read_fasta(path: Path) -> dict[str, str]:
    with open(path, "r", encoding="utf-8") as fh:
        return {r.id: str(r.seq) for r in SeqIO.parse(fh, "fasta")}


def test_fasta_record_cache_load_and_switch_files(write_fasta):
    f1 = write_fasta("one.fa", [("a", "AAAA", None), ("b", "CC", None)])
    f2 = write_fasta("two.fa", [("c", "GGG", None)])

    loc1a = recombine.RecordLocation(path=f1, description="a desc")
    loc2c = recombine.RecordLocation(path=f2, description="c desc")

    cache = recombine.FastaRecordCache()
    rec_a = cache.get("a", loc1a)
    assert rec_a.id == "a"
    assert str(rec_a.seq) == "AAAA"

    # switch file
    rec_c = cache.get("c", loc2c)
    assert rec_c.id == "c"
    assert str(rec_c.seq) == "GGG"


def test_fasta_record_cache_missing_record_raises_keyerror(write_fasta):
    f = write_fasta("one.fa", [("a", "AAAA", None)])
    loc = recombine.RecordLocation(path=f, description="x")
    cache = recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="Record 'invalid' not found"):
        cache.get("invalid", loc)


def test_fasta_record_cache_duplicate_id_within_file_raises(write_fasta):
    f = write_fasta("dup.fa", [("a", "AAAA", None), ("a", "TTTT", None)])
    cache = recombine.FastaRecordCache()
    with pytest.raises(ValueError, match="Duplicate record id 'a'"):
        cache._load_file(f)


def test_build_index_collects_locations_first_seen_and_chunks(write_fasta):
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

    locations, first_seen, chunks = recombine._build_index(CHUNK_RE, [f1, f2])

    assert "Y" in locations
    assert "X_chunk_start_0" in locations
    assert first_seen["X"] < first_seen["Y"]  # X base first observed in f1
    assert chunks["X"] == [(0, "X_chunk_start_0"), (4, "X_chunk_start_4")]


def test_build_index_duplicate_record_id_raises(write_fasta):
    f1 = write_fasta("f1.fa", [("A", "AAA", None)])
    f2 = write_fasta("f2.fa", [("A", "TTT", None)])  # duplicate id across files
    with pytest.raises(ValueError, match="Duplicate FASTA record id encountered during indexing"):
        recombine._build_index(CHUNK_RE, [f1, f2])


def test_build_index_no_headers_raises(tmp_path):
    # create an empty file that parses to no records
    p = tmp_path / "empty.fa"
    p.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="No FASTA headers found"):
        recombine._build_index(CHUNK_RE, [p])


def test_records_from_headers_passes_unchunked(write_fasta):
    f1 = write_fasta("f1.fa", [("A", "AAA", "descA")])
    locations, first_seen, chunks = recombine._build_index(CHUNK_RE, [f1])

    cache = recombine.FastaRecordCache()
    recs = list(recombine._records_from_headers(locations, first_seen, chunks, cache))
    assert [(r.id, str(r.seq)) for r in recs] == [("A", "AAA")]
    assert recs[0].description.endswith("descA")


def test_records_from_headers_reassembles_contiguous_chunks(write_fasta):
    f1 = write_fasta(
        "f1.fa",
        [
            ("X_chunk_start_0", "AAAA", "c0"),
            ("X_chunk_start_4", "TT", "c4"),
        ],
    )
    locations, first_seen, chunks = recombine._build_index(CHUNK_RE, [f1])
    cache = recombine.FastaRecordCache()

    recs = list(recombine._records_from_headers(locations, first_seen, chunks, cache))
    assert len(recs) == 1
    assert recs[0].id == "X"
    assert str(recs[0].seq) == "AAAATT"


def test_records_from_headers_noncontiguous_chunks_raises(write_fasta):
    f1 = write_fasta(
        "f1.fa",
        [
            ("X_chunk_start_0", "AAAA", None),
            ("X_chunk_start_99", "TT", None),  # wrong start
        ],
    )
    locations, first_seen, chunks = recombine._build_index(CHUNK_RE, [f1])
    cache = recombine.FastaRecordCache()

    with pytest.raises(ValueError, match="Non-contiguous chunks for 'X'"):
        list(recombine._records_from_headers(locations, first_seen, chunks, cache))


def test_records_from_headers_missing_base_record_raises():
    locations = {}
    first_seen = {"A": 0}
    chunks = {}
    cache = recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="Base record 'A' not found"):
        list(recombine._records_from_headers(locations, first_seen, chunks, cache))


def test_records_from_headers_missing_chunk_record_raises():
    # base "X" has chunks but one chunk id isn't in locations
    locations = {}
    first_seen = {"X": 0}
    chunks = {"X": [(0, "X_chunk_start_0")]}
    cache = recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="Chunk record 'X_chunk_start_0' not found"):
        list(recombine._records_from_headers(locations, first_seen, chunks, cache))


def test_agp_component_seq_plus_and_minus():
    record = SeqRecord(Seq("ACGTGGTT"), "part", "test_record")
    assert str(recombine._agp_component_seq(record, 1, 4, "+", allow_revcomp=False)) == "ACGT"
    assert str(recombine._agp_component_seq(record, 3, 6, "-", allow_revcomp=True)) == str(
        Seq("GTGG").reverse_complement()
    )

    with pytest.raises(ValueError, match="--allow-revcomp is not enabled"):
        recombine._agp_component_seq(record, 3, 6, "-", allow_revcomp=False)

    with pytest.raises(ValueError, match="Invalid AGP orientation"):
        recombine._agp_component_seq(record, 1, 2, "?", allow_revcomp=True)


def test_records_from_agp_success(write_fasta, tmp_path):
    parts = write_fasta("parts.fa", [("p1", "AAAA", "p1desc"), ("p2", "TT", "p2desc")])
    locations, _, _ = recombine._build_index(CHUNK_RE, [parts])

    # AGP describes obj = p1(1..4) + p2(1..2)
    agp = tmp_path / "x.agp"
    agp.write_text(
        "obj\t1\t4\t1\tW\tp1\t1\t4\t+\n" "obj\t5\t6\t2\tW\tp2\t1\t2\t+\n",
        encoding="utf-8",
    )
    agp_entries = parse_agp(agp, allow_revcomp=False)

    cache = recombine.FastaRecordCache()
    recs = list(recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))
    assert [(r.id, str(r.seq)) for r in recs] == [("obj", "AAAATT")]


def test_records_from_agp_sorts_parts_out_of_order(write_fasta, tmp_path):
    parts = write_fasta("parts.fa", [("p1", "AAAA", "p1desc"), ("p2", "TT", "p2desc")])
    locations, _, _ = recombine._build_index(CHUNK_RE, [parts])

    # AGP lines deliberately out-of-order (p2 line first); should still assemble p1 + p2
    agp = tmp_path / "x.agp"
    agp.write_text(
        "obj\t5\t6\t2\tW\tp2\t1\t2\t+\n" "obj\t1\t4\t1\tW\tp1\t1\t4\t+\n",
        encoding="utf-8",
    )
    agp_entries = parse_agp(agp, allow_revcomp=False)

    cache = recombine.FastaRecordCache()
    recs = list(recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))
    assert [(r.id, str(r.seq)) for r in recs] == [("obj", "AAAATT")]


def test_records_from_agp_missing_component_raises(write_fasta):
    f1 = write_fasta("parts.fa", [("p1", "AAAA", None)])
    locations, _, _ = recombine._build_index(CHUNK_RE, [f1])

    agp_entries = {
        "obj": [
            recombine.AgpEntry(
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

    cache = recombine.FastaRecordCache()
    with pytest.raises(KeyError, match="not found in indexed FASTA headers"):
        list(recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))


def test_records_from_agp_noncontiguous_object_raises(write_fasta):
    f1 = write_fasta("parts.fa", [("p1", "AAAA", None), ("p2", "TT", None)])
    locations, _, _ = recombine._build_index(CHUNK_RE, [f1])

    agp_entries = {
        "obj": [
            recombine.AgpEntry("obj", 1, 4, 1, "p1", 1, 4, "+"),
            # should start at 5, but starts at 6
            recombine.AgpEntry("obj", 6, 7, 2, "p2", 1, 2, "+"),
        ]
    }
    cache = recombine.FastaRecordCache()
    with pytest.raises(ValueError, match="Non-contiguous AGP"):
        list(recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))


def test_records_from_agp_length_mismatch_raises(write_fasta):
    f1 = write_fasta("parts.fa", [("p1", "AAAA", None)])
    locations, _, _ = recombine._build_index(CHUNK_RE, [f1])

    # asks for 1..10 but record shorter => extracted length != expected_len
    agp_entries = {"obj": [recombine.AgpEntry("obj", 1, 10, 1, "p1", 1, 10, "+")]}
    cache = recombine.FastaRecordCache()
    with pytest.raises(ValueError, match="Length mismatch"):
        list(recombine._records_from_agp(agp_entries, locations, cache, allow_revcomp=False))


def test_recombine_fasta_header_mode_end_to_end_manifest_order(write_fasta, tmp_path):
    f1 = write_fasta("f1.fa", [("Y", "CC", None)])
    f2 = write_fasta("f2.fa", [("X_chunk_start_0", "AAAA", None), ("X_chunk_start_4", "TT", None)])

    manifest = write_manifest(tmp_path / "manifest.txt", [str(f1), str(f2)])

    out = tmp_path / "out.fa"
    recombine.recombine_fasta(
        fasta_manifest=manifest,
        out_fasta=out,
        chunk_re=CHUNK_RE,
        agp_file=None,
        allow_revcomp=False,
    )

    seqs = read_fasta(out)
    assert seqs["X"] == "AAAATT"
    assert seqs["Y"] == "CC"

    ids_in_file = [r.id for r in SeqIO.parse(out, "fasta")]
    assert ids_in_file == ["Y", "X"]


def test_recombine_fasta_agp_mode_end_to_end(write_fasta, tmp_path):
    parts = write_fasta("parts.fa", [("p1", "AACCGG", None), ("p2", "TTAA", None)])

    agp = tmp_path / "x.agp"
    agp.write_text(
        "obj\t1\t2\t1\tW\tp1\t1\t2\t+\n" "obj\t3\t4\t2\tW\tp2\t3\t4\t-\n",  # part2 slice 'AA' revcomp => 'TT'
        encoding="utf-8",
    )

    manifest = write_manifest(tmp_path / "manifest.txt", [str(parts)])

    out = tmp_path / "out.fa"
    recombine.recombine_fasta(
        fasta_manifest=manifest,
        out_fasta=out,
        chunk_re=CHUNK_RE,
        agp_file=agp,
        allow_revcomp=True,
    )
    seqs = read_fasta(out)
    assert seqs["obj"] == "AATT"


def test_recombine_header_mode_with_alternative_regex_end_to_end(write_fasta, tmp_path: Path):
    alt_chunk_re = validate_regex(r"^(?P<base>.+)\.chunk\.(?P<start>\d+)$")

    f1 = write_fasta("d/a.fa", [("X.chunk.0", "AAAA", None)])
    f2 = write_fasta("d/b.fa", [("X.chunk.4", "TT", None)])
    f3 = write_fasta("d/c.fa", [("Y", "CC", None)])

    manifest = write_manifest(tmp_path / "manifest.txt", [str(f1), str(f2), str(f3)])

    out = tmp_path / "out.fa"
    recombine.recombine_fasta(
        fasta_manifest=manifest,
        out_fasta=out,
        chunk_re=alt_chunk_re,
        agp_file=None,
        allow_revcomp=False,
    )

    seqs = read_fasta(out)
    assert seqs["X"] == "AAAATT"
    assert seqs["Y"] == "CC"


def test_parse_args_minimal_required(write_fasta, tmp_path):
    f = write_fasta("a.fa", [("A", "AAA", None)])
    manifest = write_manifest(tmp_path / "manifest.txt", [str(f)])
    out = tmp_path / "out.fa"
    args = recombine.parse_args(
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
    assert validate_regex(args.chunk_id_regex).match("X_chunk_start_0")


def test_parse_args_boolean_flags_manifest(write_fasta, tmp_path):
    f = write_fasta("a.fa", [("A", "AAA", None)])
    manifest = write_manifest(tmp_path / "manifest.txt", [str(f)])
    out = tmp_path / "out.fa"

    args = recombine.parse_args(
        [
            "--fasta-manifest",
            str(manifest),
            "--out-fasta",
            str(out),
            "--allow-revcomp",
        ]
    )
    assert args.allow_revcomp is True
