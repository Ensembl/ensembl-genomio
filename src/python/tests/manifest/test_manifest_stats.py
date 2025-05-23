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
"""Unit testing of `ensembl.io.genomio.manifest.manifest_stats` module."""

from contextlib import nullcontext as no_raise
from pathlib import Path
from shutil import copy
from typing import ContextManager

import pytest

from ensembl.io.genomio.manifest import Manifest
from ensembl.io.genomio.manifest.manifest_stats import ManifestStats
from ensembl.io.genomio.utils import print_json


@pytest.fixture(name="manifest_path")
def fixture_manifest_path(data_dir: Path) -> Path:
    """Manifest dir and files with all expected cases."""
    return data_dir / "full_data" / "manifest.json"


def test_manifest_stats_init(manifest_path: Path) -> None:
    """Test `ManifestStats.__init__()`."""
    manifest_stats = ManifestStats(manifest_path)
    assert manifest_stats


def test_add_error(manifest_path: Path) -> None:
    """Test `ManifestStats.has_lengths()`."""
    stats = ManifestStats(manifest_path)
    assert not stats.errors
    stats.add_error("lorem")
    assert stats.errors[0] == "lorem"


@pytest.mark.parametrize(
    ("manifest_dir", "expected_lengths", "expected_circular"),
    [
        pytest.param("", {}, {}),
        pytest.param(
            "seq_regions",
            {"seq1": 100, "seq1_gb": 100, "seq1_insdc": 100, "seq2": 10},
            {"seq1": True, "seq2": False},
        ),
    ],
)
def test_load_seq_regions(
    data_dir: Path,
    tmp_path: Path,
    manifest_dir: str,
    expected_lengths: dict[str, int],
    expected_circular: dict[str, bool],
) -> None:
    """Test `ManifestStats.load_seq_regions()`.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Fixture which will provide a temporary directory unique to the test invocation.
        manifest_dir: Directory name with the manifest file in it.
        expected_lengths: Expected length data from the files.
        expected_circular: Expected circular data from the files.

    """
    if manifest_dir:
        seq_manifest = data_dir / manifest_dir / "manifest.json"
    else:
        seq_manifest = tmp_path / "manifest.json"
        with seq_manifest.open("w") as manifest_fh:
            manifest_fh.write("{}")

    stats = ManifestStats(data_dir / seq_manifest)
    stats.load_seq_regions()
    assert stats.lengths["seq_regions"] == expected_lengths
    assert stats.circular["seq_regions"] == expected_circular


@pytest.mark.parametrize(
    ("fasta_str", "expected_lengths", "expected_error"),
    [
        pytest.param(">prot1\nCGTA\n", {"prot1": 4}, "", id="Normal DNA seq"),
        pytest.param("", {}, "", id="Empty fasta"),
        pytest.param("CGTA", {}, "No sequences found", id="No sequences in fasta"),
        pytest.param(">\nCGTA\n", {}, "1 sequences with empty ids", id="empty ID"),
    ],
)
def test_load_dna_fasta_lengths(
    tmp_path: Path,
    fasta_str: str,
    expected_lengths: dict[str, int],
    expected_error: str,
) -> None:
    """Test `ManifestStats.load_dna_fasta_lengths()`.

    Args:
        tmp_path: Fixture which will provide a temporary directory unique to the test invocation.
        fasta_str: Content of a test input FASTA DNA file.
        expected_lengths: Expected length data from the files.
        expected_error: Expected errors while loading.

    """
    fasta_path = tmp_path / "fasta_dna.fasta"
    with fasta_path.open("w") as fasta_fh:
        fasta_fh.write(fasta_str)
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.root_dir / "manifest.json")
    stats.load_dna_fasta_lengths()
    assert stats.lengths["dna_sequences"] == expected_lengths
    if expected_error == "":
        assert not stats.errors
    else:
        assert len(stats.errors) == 1
        assert stats.errors[0].startswith(expected_error)


@pytest.mark.parametrize(
    ("fasta_str", "ignore_final_stops", "expected_lengths", "expected_error"),
    [
        pytest.param(">prot1\nMAGIC\n", False, {"prot1": 5}, "", id="Normal prot"),
        pytest.param("", False, {}, "", id="Empty fasta"),
        pytest.param("AH", False, {}, "No sequences found", id="No sequences in fasta"),
        pytest.param(">\nMAGIC\n", False, {}, "1 sequences with empty ids", id="Empty ID"),
        pytest.param(
            ">prot1\nMAGIC*\n", False, {"prot1": 6}, "1 sequences with stop codons", id="End stop codon"
        ),
        pytest.param(">prot1\nMAGIC\n", True, {"prot1": 5}, "", id="End stop codon, ignore"),
        pytest.param(
            ">prot1\nMAG*IC\n", False, {"prot1": 6}, "1 sequences with stop codons", id="In stop codon"
        ),
        pytest.param(
            ">prot1\nMAGIC\n>prot1\nGICMA\n",
            False,
            {"prot1": 5},
            "1 non unique sequence ids",
            id="Duplicate ID",
        ),
    ],
)
def test_load_peptides_fasta_lengths(
    tmp_path: Path,
    fasta_str: str,
    ignore_final_stops: bool,
    expected_lengths: dict[str, int],
    expected_error: str,
) -> None:
    """Test `ManifestStats.load_peptides_fasta_lengths()`.

    Args:
        tmp_path: Fixture which will provide a temporary directory unique to the test invocation.
        fasta_str: Content of a test input FASTA DNA file.
        ignore_final_stops: Ignore final stops for the protein sequences.
        expected_lengths: Expected length data from the files.
        expected_error: Expected errors while loading.

    """
    fasta_path = tmp_path / "fasta_pep.fasta"
    with fasta_path.open("w") as fasta_fh:
        fasta_fh.write(fasta_str)
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.root_dir / "manifest.json", ignore_final_stops=ignore_final_stops)
    stats.load_peptides_fasta_lengths()
    assert stats.lengths["peptide_sequences"] == expected_lengths
    if expected_error == "":
        assert not stats.errors
    else:
        assert len(stats.errors) == 1
        assert stats.errors[0].startswith(expected_error)


@pytest.mark.parametrize(
    ("json_data", "expected_key", "expected_data"),
    [
        pytest.param(None, "", {}, id="No JSON"),
        pytest.param({}, "", {}, id="Empty JSON"),
        pytest.param({"id": "gene1", "object_type": "gene"}, "ann_genes", {"gene1": 1}, id="1 gene"),
        pytest.param(
            {"id": "pep1", "object_type": "translation"},
            "ann_translations",
            {"pep1": 1},
            id="1 translation",
        ),
        pytest.param(
            {"id": "te1", "object_type": "transposable_element"},
            "ann_transposable_elements",
            {"te1": 1},
            id="1 transposable element",
        ),
    ],
)
def test_load_functional_annotation(
    tmp_path: Path,
    json_data: dict | None,
    expected_key: str,
    expected_data: dict[str, int],
) -> None:
    """Test `ManifestStats.load_functional_annotation()`.

    Args:
        tmp_path: Fixture which will provide a temporary directory unique to the test invocation.
        json_data: JSON data to write to the functional annotation file.
        expected_key: Key to check in the loaded functional annotation.
        expected_data: Expected length data for the functional annotation key.

    """
    func_path = tmp_path / "functional_annotation.json"
    if json_data is not None:
        if json_data:
            print_json(func_path, [json_data])
        else:
            print_json(func_path, [])
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.root_dir / "manifest.json")
    stats.load_functional_annotation()
    if expected_key:
        assert stats.lengths[expected_key] == expected_data


@pytest.mark.parametrize(
    ("gff3_path", "expected_data"),
    [
        pytest.param(None, {}, id="No GFF3"),
        pytest.param(
            "get_gff3/ok_gene.gff",
            {
                "gff3_seq_regions": {"scaffold1": 1000},
                "gff3_genes": {"LOREMIPSUM1": 1000},
                "gff3_all_translations": {"LOREMIPSUM1_t1_cds": 266},
                "gff3_translations": {"LOREMIPSUM1_t1_cds": 266},
            },
            id="1 gene, 1 translation",
        ),
        pytest.param(
            "get_gff3/te_gene.gff",
            {
                "gff3_transposable_elements": {"LOREMIPSUM1": 1000},
            },
            id="1 TE",
        ),
        pytest.param(
            "get_gff3/pseudogene_cds.gff",
            {
                "gff3_genes": {"LOREMIPSUM1": 1000},
                "gff3_translations": {},
                "gff3_all_translations": {"LOREMIPSUM1_t1_cds": 266},
            },
            id="1 pseudogene with CDS",
        ),
    ],
)
def test_load_gff3(
    tmp_path: Path,
    data_dir: Path,
    gff3_path: str | None,
    expected_data: dict[str, int],
) -> None:
    """Test `ManifestStats.load_gff3()`.

    Args:
        tmp_path: Fixture which will provide a temporary directory unique to the test invocation.
        data_dir: Module's test data directory fixture.
        gff3_path: Path to a test GFF3 file.
        expected_data: Expected length data extracted from the GFF3 file.

    """
    if gff3_path:
        copy(data_dir / gff3_path, tmp_path / "test.gff3")
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.root_dir / "manifest.json")
    stats.load_gff3()
    expected_lengths = {**stats.lengths, **expected_data}
    assert stats.lengths == expected_lengths


def test_load_genome(tmp_path: Path) -> None:
    """Test `ManifestStats.load_genome()`."""
    genome_data = {"LOREM": "1"}
    genome_path = tmp_path / "genome.json"
    print_json(genome_path, genome_data)
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.file_path)
    stats.load_genome()
    assert stats.genome == genome_data


def test_prepare_integrity_data(tmp_path: Path) -> None:
    """Test `ManifestStats.prepare_integrity_data()`."""
    manifest_path = tmp_path / "manifest.json"
    print_json(manifest_path, {})
    stats = ManifestStats(manifest_path)
    stats.prepare_integrity_data()
    assert stats.lengths == {
        "agp": {},
        "ann_genes": {},
        "ann_translations": {},
        "ann_transposable_elements": {},
        "annotations": {},
        "dna_sequences": {},
        "gff3_all_translations": {},
        "gff3_genes": {},
        "gff3_seq_regions": {},
        "gff3_translations": {},
        "gff3_transposable_elements": {},
        "peptide_sequences": {},
        "seq_region_levels": {},
        "seq_regions": {},
    }


@pytest.mark.parametrize(
    ("key", "expected_data", "expectation"),
    [
        pytest.param("ann_genes", True, no_raise(), id="Loaded func annotation"),
        pytest.param("gff3_genes", False, no_raise(), id="Not loaded gff genes"),
        pytest.param("foobar", False, pytest.raises(KeyError), id="Unknown key"),
    ],
)
def test_has_lengths(
    manifest_path: Path,
    key: str,
    expected_data: dict[str, bool],
    expectation: ContextManager,
) -> None:
    """Test `ManifestStats.has_lengths()`.

    Args:
        manifest_path: Path to the manifest file.
        key: Key for the length dict.
        expected_data: If lengths exist for that key.
        expectation: Expected exception.

    """
    stats = ManifestStats(manifest_path)
    stats.load_functional_annotation()
    with expectation:
        assert stats.has_lengths(key) == expected_data


@pytest.mark.parametrize(
    ("key", "expected_data", "expectation"),
    [
        pytest.param("ann_genes", {"ECC02_000372": 1}, no_raise(), id="OK"),
        pytest.param("foobar", {}, pytest.raises(KeyError), id="Unknown key"),
    ],
)
def test_get_lengths(
    manifest_path: Path,
    key: str,
    expected_data: dict[str, bool],
    expectation: ContextManager,
) -> None:
    """Test `ManifestStats.get_lengths()`.

    Args:
        manifest_path: Path to the manifest file.
        key: Key for the length dict.
        expected_data: Expected length information for that key.
        expectation: Expected exception.

    """
    stats = ManifestStats(manifest_path)
    stats.load_functional_annotation()
    with expectation:
        assert stats.get_lengths(key) == expected_data


@pytest.mark.parametrize(
    ("key", "expected_data", "expectation"),
    [
        pytest.param(
            "seq_regions", {"JABDHM010000001.1": True, "JABDHM010000002.1": False}, no_raise(), id="OK"
        ),
        pytest.param("foobar", {}, pytest.raises(KeyError), id="Unknown key"),
    ],
)
def test_get_circular(
    manifest_path: Path,
    key: str,
    expected_data: dict[str, bool],
    expectation: ContextManager,
) -> None:
    """Test `ManifestStats.get_circular()`.

    Args:
        manifest_path: Path to the manifest file.
        key: Key for the circular dict.
        expected_data: Expected circular information for that key.
        expectation: Expected exception.

    """
    stats = ManifestStats(manifest_path)
    stats.load_seq_regions()

    with expectation:
        assert stats.get_circular(key) == expected_data
