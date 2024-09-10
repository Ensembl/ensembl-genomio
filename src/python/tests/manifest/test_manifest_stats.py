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

from pathlib import Path

import pytest
from pytest import raises, param, mark

from ensembl.io.genomio.manifest import Manifest
from ensembl.io.genomio.manifest.manifest_stats import ManifestStats
from ensembl.io.genomio.utils import print_json


@pytest.fixture(name="manifest_path")
def fixture_manifest_path(data_dir: Path) -> Path:
    """Manifest dir and files with all expected cases."""
    return data_dir / "full_data/manifest.json"


def test_manifest_stats_init(manifest_path: Path) -> None:
    """Tests `ManifestStats.__init__()`."""
    manifest_stats = ManifestStats(manifest_path)
    assert manifest_stats


def test_add_error(manifest_path: Path) -> None:
    """Tests `ManifestStats.has_lengths()`."""
    stats = ManifestStats(manifest_path)
    assert not stats.errors
    stats.add_error("lorem")
    assert stats.errors[0] == "lorem"


@mark.parametrize(
    "manifest_dir, expected_lengths, expected_circular",
    [
        param("", {}, {}),
        param(
            "seq_regions",
            {"seq1": 100, "seq1_gb": 100, "seq1_insdc": 100, "seq2": 10},
            {"seq1": True, "seq2": False},
        ),
    ],
)
def test_get_seq_regions(
    data_dir: Path,
    tmp_path: Path,
    manifest_dir: str,
    expected_lengths: dict[str, int],
    expected_circular: dict[str, bool],
) -> None:
    """Tests `ManifestStats.get_seq_regions()`."""
    if manifest_dir:
        seq_manifest = data_dir / manifest_dir / "manifest.json"
    else:
        seq_manifest = tmp_path / "manifest.json"
        with seq_manifest.open("w") as manifest_fh:
            manifest_fh.write("{}")

    stats = ManifestStats(data_dir / seq_manifest)
    stats.get_seq_regions()
    assert stats.lengths["seq_regions"] == expected_lengths
    assert stats.circular["seq_regions"] == expected_circular


@mark.parametrize(
    "fasta_str, expected_lengths, expected_error",
    [
        param(">prot1\nCGTA\n", {"prot1": 4}, "", id="Normal DNA seq"),
        param("", {}, "", id="Empty fasta"),
        param("CGTA", {}, "No sequences found", id="No sequences in fasta"),
        param(">\nCGTA\n", {}, "1 sequences with empty ids", id="empty ID"),
    ],
)
def test_get_dna_fasta_lengths(
    tmp_path: Path,
    fasta_str: str,
    expected_lengths: dict[str, int],
    expected_error: str,
) -> None:
    """Tests `ManifestStats.get_dna_fasta_lengths()`."""
    fasta_path = tmp_path / "fasta_dna.fasta"
    with fasta_path.open("w") as fasta_fh:
        fasta_fh.write(fasta_str)
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.dir / "manifest.json")
    stats.get_dna_fasta_lengths()
    assert stats.lengths["dna_sequences"] == expected_lengths
    if expected_error == "":
        assert not stats.errors
    else:
        assert len(stats.errors) == 1
        assert stats.errors[0].startswith(expected_error)


@mark.parametrize(
    "fasta_str, ignore_final_stops, expected_lengths, expected_error",
    [
        param(">prot1\nMAGIC\n", False, {"prot1": 5}, "", id="Normal prot"),
        param("", False, {}, "", id="Empty fasta"),
        param("AH", False, {}, "No sequences found", id="No sequences in fasta"),
        param(">\nMAGIC\n", False, {}, "1 sequences with empty ids", id="empty ID"),
        param(">prot1\nMAGIC*\n", False, {"prot1": 6}, "1 sequences with stop codons", id="End stop codon"),
        param(">prot1\nMAGIC\n", True, {"prot1": 5}, "", id="End stop codon, ignore"),
        param(">prot1\nMAG*IC\n", False, {"prot1": 6}, "1 sequences with stop codons", id="In stop codon"),
        param(
            ">prot1\nMAGIC\n>prot1\nGICMA\n",
            False,
            {"prot1": 5},
            "1 non unique sequence ids",
            id="duplicate ID",
        ),
    ],
)
def test_get_peptides_fasta_lengths(
    tmp_path: Path,
    fasta_str: str,
    ignore_final_stops: bool,
    expected_lengths: dict[str, int],
    expected_error: str,
) -> None:
    """Tests `ManifestStats.get_peptides_fasta_lengths()`."""
    fasta_path = tmp_path / "fasta_pep.fasta"
    with fasta_path.open("w") as fasta_fh:
        fasta_fh.write(fasta_str)
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.dir / "manifest.json")
    stats.ignore_final_stops = ignore_final_stops
    stats.get_peptides_fasta_lengths()
    assert stats.lengths["peptide_sequences"] == expected_lengths
    if expected_error == "":
        assert not stats.errors
    else:
        assert len(stats.errors) == 1
        assert stats.errors[0].startswith(expected_error)


@mark.parametrize(
    "json_data, expected_key, expected_data",
    [
        param(None, "", {}, id="No JSON"),
        param({}, "", {}, id="Empty JSON"),
        param({"id": "gene1", "object_type": "gene"}, "ann_genes", {"gene1": 1}, id="1 gene"),
        param(
            {"id": "pep1", "object_type": "translation"}, "ann_translations", {"pep1": 1}, id="1 translation"
        ),
        param(
            {"id": "te1", "object_type": "transposable_element"},
            "ann_transposable_elements",
            {"te1": 1},
            id="1 te",
        ),
    ],
)
def test_get_functional_annotations(
    tmp_path: Path,
    json_data: dict | None,
    expected_key: str,
    expected_data: dict[str, int],
) -> None:
    """Tests `ManifestStats.get_functional_annotation()`."""
    func_path = tmp_path / "functional_annotation.json"
    if json_data is not None:
        if json_data:
            print_json(func_path, [json_data])
        else:
            print_json(func_path, [])
    manifest = Manifest(tmp_path)
    manifest.create()
    stats = ManifestStats(manifest.dir / "manifest.json")
    stats.get_functional_annotation()
    if expected_key:
        assert stats.lengths[expected_key] == expected_data


def test_has_lengths(manifest_path: Path) -> None:
    """Tests `ManifestStats.has_lengths()`."""
    stats = ManifestStats(manifest_path)
    with raises(KeyError):
        stats.has_lengths("foobar")
    assert not stats.has_lengths("ann_genes")
    stats.get_functional_annotation()
    assert stats.has_lengths("ann_genes")


def test_get_lengths(manifest_path: Path) -> None:
    """Tests `ManifestStats.get_lengths()`."""
    stats = ManifestStats(manifest_path)
    with raises(KeyError):
        stats.get_lengths("foobar")
    assert not stats.has_lengths("ann_genes")
    stats.get_functional_annotation()
    circular = stats.get_lengths("ann_genes")
    assert circular == {"ECC02_000372": 1}


def test_get_circular(manifest_path: Path) -> None:
    """Tests `ManifestStats.get_circular()`."""
    stats = ManifestStats(manifest_path)
    with raises(KeyError):
        stats.get_circular("foobar")
    assert not stats.get_circular("seq_regions")
    stats.get_seq_regions()
    lengths = stats.get_circular("seq_regions")
    assert lengths == {"JABDHM010000001.1": True, "JABDHM010000002.1": False}
