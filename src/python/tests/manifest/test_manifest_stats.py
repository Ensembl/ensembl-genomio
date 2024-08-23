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
from pytest import raises

from ensembl.io.genomio.manifest.manifest_stats import ManifestStats, InvalidIntegrityError


def test_manifest_stats_init(data_dir: Path) -> None:
    """Tests `ManifestStats.__init__()`."""
    manifest_stats = ManifestStats(data_dir / "full_data/manifest.json")
    assert manifest_stats


def test_has_lengths(data_dir: Path):
    """Tests `ManifestStats.has_lengths()`."""
    stats = ManifestStats(data_dir / "full_data/manifest.json")
    with raises(KeyError):
        stats.has_lengths("foobar")
    assert not stats.has_lengths("ann_genes")
    stats.get_functional_annotation(data_dir / "full_data" / "functional_annotation.json")
    assert stats.has_lengths("ann_genes")


def test_get_lengths(data_dir: Path):
    """Tests `ManifestStats.get_lengths()`."""
    stats = ManifestStats(data_dir / "full_data/manifest.json")
    with raises(KeyError):
        stats.get_lengths("foobar")
    assert not stats.has_lengths("ann_genes")
    stats.get_functional_annotation(data_dir / "full_data" / "functional_annotation.json")
    circular = stats.get_lengths("ann_genes")
    assert circular == {"ECC02_000372": 1}


def test_get_circular(data_dir: Path):
    """Tests `ManifestStats.get_circular()`."""
    stats = ManifestStats(data_dir / "full_data/manifest.json")
    with raises(KeyError):
        stats.get_circular("foobar")
    assert not stats.get_circular("seq_regions")
    stats.get_seq_regions()
    lengths = stats.get_circular("seq_regions")
    assert lengths == {"JABDHM010000001.1": True, "JABDHM010000002.1": False}
