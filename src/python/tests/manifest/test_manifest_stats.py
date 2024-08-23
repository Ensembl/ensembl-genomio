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

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager

import pytest
from pytest import raises, param, mark

from ensembl.io.genomio.manifest.manifest_stats import ManifestStats, InvalidIntegrityError


@pytest.fixture(name="manifest_path")
def fixture_manifest_path(data_dir: Path) -> Path:
    """Manifest dir and files with all expected cases."""
    return data_dir / "full_data/manifest.json"


def test_manifest_stats_init(manifest_path: Path) -> None:
    """Tests `ManifestStats.__init__()`."""
    manifest_stats = ManifestStats(manifest_path)
    assert manifest_stats


def test_add_error(manifest_path: Path):
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
):
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


def test_has_lengths(manifest_path: Path):
    """Tests `ManifestStats.has_lengths()`."""
    stats = ManifestStats(manifest_path)
    with raises(KeyError):
        stats.has_lengths("foobar")
    assert not stats.has_lengths("ann_genes")
    stats.get_functional_annotation(manifest_path.parent / "functional_annotation.json")
    assert stats.has_lengths("ann_genes")


def test_get_lengths(manifest_path: Path):
    """Tests `ManifestStats.get_lengths()`."""
    stats = ManifestStats(manifest_path)
    with raises(KeyError):
        stats.get_lengths("foobar")
    assert not stats.has_lengths("ann_genes")
    stats.get_functional_annotation(manifest_path.parent / "functional_annotation.json")
    circular = stats.get_lengths("ann_genes")
    assert circular == {"ECC02_000372": 1}


def test_get_circular(manifest_path: Path):
    """Tests `ManifestStats.get_circular()`."""
    stats = ManifestStats(manifest_path)
    with raises(KeyError):
        stats.get_circular("foobar")
    assert not stats.get_circular("seq_regions")
    stats.get_seq_regions()
    lengths = stats.get_circular("seq_regions")
    assert lengths == {"JABDHM010000001.1": True, "JABDHM010000002.1": False}
