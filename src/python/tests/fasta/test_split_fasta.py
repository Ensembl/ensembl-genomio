# tests/test_split_fasta.py
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.io.genomio.fasta.split_fasta import main as SplitFasta


def write_fasta(path: Path, records):
    """Write a list of SeqRecord objects to a FASTA file."""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        SeqIO.write(records, fh, "fasta")


def list_output_fastas(out_dir: Path):
    """Return all FASTA files produced under the output directory."""
    return sorted(out_dir.rglob("*.fa"))


def read_all_ids_from_fastas(out_dir: Path):
    """Read and return all sequence IDs from all FASTA files under out_dir."""
    ids = []
    for fa in list_output_fastas(out_dir):
        with open(fa, "r", encoding="utf-8") as fh:
            ids.extend([r.id for r in SeqIO.parse(fh, "fasta")])
    return ids


def parse_agp_lines(agp_path: Path):
    """
    Parse an AGP file into a list of column lists, excluding comments
    and blank lines.
    """
    lines = [l.rstrip("\n") for l in agp_path.read_text(encoding="utf-8").splitlines()]
    lines = [l for l in lines if l and not l.startswith("#")]
    return [l.split("\t") for l in lines]


def test_no_agp_by_default(tmp_path: Path):
    """
    By default, splitting a FASTA should produce one or more FASTA outputs
    but must NOT create an AGP file unless write_agp is explicitly enabled.
    """
    input_fasta = tmp_path / "in.fa"
    out = tmp_path / "out"
    write_fasta(input_fasta, [SeqRecord(Seq("ACGT"), id="seq1", description="")])

    SplitFasta(
        [
            "--fasta-file",
            str(input_fasta),
            "--out-dir",
            str(out),
        ]
    )

    assert not (out / "in.agp").exists()
    assert len(list_output_fastas(out)) >= 1


def test_split_by_max_seqs_per_file(tmp_path: Path):
    """
    When max_seqs_per_file is set, sequences should be split across
    multiple FASTA files while preserving original sequence order
    and IDs.
    """
    input_fasta = tmp_path / "in.fa"
    out = tmp_path / "out"
    recs = [
        SeqRecord(Seq("A" * 10), id="s1", description=""),
        SeqRecord(Seq("C" * 10), id="s2", description=""),
        SeqRecord(Seq("G" * 10), id="s3", description=""),
    ]
    write_fasta(input_fasta, recs)

    SplitFasta(
        [
            "--fasta-file",
            str(input_fasta),
            "--out-dir",
            str(out),
            "--max-seqs-per-file",
            "2",
            "--write-agp",
        ]
    )
    fas = list_output_fastas(out)
    assert len(fas) == 2
    assert read_all_ids_from_fastas(out) == ["s1", "s2", "s3"]


def test_chunk_merge_final_small_chunk_and_agp(tmp_path: Path):
    """
    When force_max_seq_length is enabled, long sequences are chunked.
    If the final chunk is shorter than min_chunk_length, it should be
    merged with the previous chunk, and the AGP file must reflect the
    merged coordinates correctly.
    """
    input_fasta = tmp_path / "in.fa"
    out = tmp_path / "out"
    write_fasta(input_fasta, [SeqRecord(Seq("A" * 2100), id="chr1", description="chr1")])

    SplitFasta(
        [
            "--fasta-file",
            str(input_fasta),
            "--out-dir",
            str(out),
            "--write-agp",
            "--force-max-seq-length",
            "--max-seq-length-per-file",
            "1000",
            "--min-chunk-length",
            "200",
            "--max-seqs-per-file",
            "100000",  # avoid seq-count splitting interfering
        ]
    )

    # 2 chunks expected after merge
    assert read_all_ids_from_fastas(out) == [
        "chr1_chunk_start_0",
        "chr1_chunk_start_1000",
    ]

    agp = out / "in.agp"
    assert agp.exists()

    cols = parse_agp_lines(agp)
    assert len(cols) == 2

    # object, obj_start, obj_end, part_no, type, comp_id, comp_start, comp_end, orientation
    assert cols[0] == [
        "chr1",
        "1",
        "1000",
        "1",
        "W",
        "chr1_chunk_start_0",
        "1",
        "1000",
        "+",
    ]
    assert cols[1] == [
        "chr1",
        "1001",
        "2100",
        "2",
        "W",
        "chr1_chunk_start_1000",
        "1",
        "1100",
        "+",
    ]


def test_agp_part_numbers_restart_per_object(tmp_path: Path):
    """
    AGP part numbers must restart at 1 for each new input sequence
    (object), even when multiple sequences are chunked in the same run.
    """
    input_fasta = tmp_path / "in.fa"
    out = tmp_path / "out"
    recs = [
        SeqRecord(Seq("A" * 1200), id="obj1", description=""),
        SeqRecord(Seq("C" * 1200), id="obj2", description=""),
    ]
    write_fasta(input_fasta, recs)

    SplitFasta(
        [
            "--fasta-file",
            str(input_fasta),
            "--out-dir",
            str(out),
            "--write-agp",
            "--force-max-seq-length",
            "--max-seq-length-per-file",
            "1000",
            "--min-chunk-length",
            "100",  # => 2 chunks each, no merge
        ]
    )

    cols = parse_agp_lines(out / "in.agp")

    by_obj = {}
    for c in cols:
        by_obj.setdefault(c[0], []).append(int(c[3]))

    assert by_obj["obj1"] == [1, 2]
    assert by_obj["obj2"] == [1, 2]
