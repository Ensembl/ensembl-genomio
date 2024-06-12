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
"""Unit testing of `ensembl.io.genomio.genome_stats.compare` module.

Typical usage example::
    $ pytest test_compare.py

"""

from pathlib import Path
from typing import Callable, Dict

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_stats import compare


@pytest.mark.dependency(name="test_compare_dicts")
@pytest.mark.parametrize(
    "ncbi, core, output",
    [
        pytest.param({}, {}, {}, id="empty_dicts"),
        pytest.param({"a": 0}, {"a": 0}, {}, id="same_dicts_zero_values"),
        pytest.param({"a": 3}, {"a": 3}, {"same": {"a": 3}}, id="same_dicts_non_zero_values"),
        pytest.param(
            {"a": 3}, {"a": 5}, {"different": {"a": {"ncbi": 3, "core": 5, "diff": 2}}}, id="different_dicts"
        ),
        pytest.param(
            {"a": 3, "b": 5},
            {"a": 3, "b": 3},
            {"same": {"a": 3}, "different": {"b": {"ncbi": 5, "core": 3, "diff": -2}}},
            id="partially_similar_dicts",
        ),
    ],
)
def test_compare_dicts(ncbi: Dict[str, int], core: Dict[str, int], output: Dict[str, Dict]) -> None:
    """Tests the `compare._compare_dicts()` method.

    Args:
        ncbi: NCBI dataset statistics in key-value pairs.
        core: Core database statistics in key-value pairs.
        output: Expected output when comparing both dictionaries.

    """
    result = compare.stats_dict_cmp(ncbi, core)
    assert not DeepDiff(result, output)


@pytest.mark.dependency(name="test_compare_assembly", depends=["test_compare_dicts"])
@pytest.mark.parametrize(
    "ncbi_file, core_file, output_file",
    [
        ("ncbi_unannotated.json", "core_unannotated.json", "output_unannotated.json"),
        ("ncbi_annotated.json", "core_annotated.json", "output_annotated.json"),
    ],
)
def test_compare_assembly(json_data: Callable, ncbi_file: Dict, core_file: Dict, output_file: Dict) -> None:
    """Tests the `compare.compare_assembly()` method.

    Args:
        json_data: JSON test file parsing fixture
        ncbi_file: NCBI dataset assembly statistics file.
        core_file: Core database assembly statistics file.
        output_file: Expected output when comparing both sources.

    """
    ncbi_stats = json_data(ncbi_file)["reports"][0]
    core_stats = json_data(core_file)["assembly_stats"]
    output = json_data(output_file)["assembly_diff"]
    result = compare.compare_assembly(ncbi_stats, core_stats)
    assert not DeepDiff(result, output)


@pytest.mark.dependency(name="test_compare_annotation", depends=["test_compare_dicts"])
@pytest.mark.parametrize(
    "ncbi_file, core_file, output_file",
    [
        ("ncbi_unannotated.json", "core_unannotated.json", "output_unannotated.json"),
        ("ncbi_annotated.json", "core_annotated.json", "output_annotated.json"),
    ],
)
def test_compare_annotation(json_data: Callable, ncbi_file: str, core_file: str, output_file: str) -> None:
    """Tests the `compare.compare_annotation()` method.

    Args:
        json_data: JSON test file parsing fixture
        ncbi_file: NCBI dataset annotation statistics file.
        core_file: Core database annotation statistics file.
        output_file: Expected output when comparing both sources.

    """
    ncbi_stats = json_data(ncbi_file)["reports"][0].get("annotation_info", {})
    if ncbi_stats:
        ncbi_stats = ncbi_stats["stats"]["gene_counts"]
    core_stats = json_data(core_file).get("annotation_stats", {})
    output = json_data(output_file).get("annotation_diff", {})
    result = compare.compare_annotation(ncbi_stats, core_stats)
    assert not DeepDiff(result, output)


@pytest.mark.dependency(
    name="test_compare_stats", depends=["test_compare_assembly", "test_compare_annotation"]
)
@pytest.mark.parametrize(
    "ncbi_file, core_file, output_file",
    [
        ("ncbi_unannotated.json", "core_unannotated.json", "output_unannotated.json"),
        ("ncbi_annotated.json", "core_annotated.json", "output_annotated.json"),
    ],
)
def test_compare_stats(json_data: Callable, ncbi_file: str, core_file: str, output_file: str) -> None:
    """Tests the `compare.compare_stats()` method.

    Args:
        json_data: JSON test file parsing fixture
        ncbi_file: NCBI dataset genome statistics file.
        core_file: Core database genome statistics file.
        output_file: Expected output when comparing both sources.

    """
    ncbi_stats = json_data(ncbi_file)["reports"][0]
    core_stats = json_data(core_file)
    output = json_data(output_file)
    result = compare.compare_stats(ncbi_stats, core_stats)
    assert not DeepDiff(result, output)


@pytest.mark.dependency(name="test_compare_stats_files", depends=["test_compare_stats"])
@pytest.mark.parametrize(
    "ncbi_file, core_file, output_file",
    [
        ("ncbi_annotated.json", "core_annotated.json", "output_annotated.json"),
    ],
)
def test_compare_stats_files(
    data_dir: Path, json_data: Callable, ncbi_file: str, core_file: str, output_file: str
) -> None:
    """Tests the `compare.compare_stats_files()` method.

    Args:
        data_dir: Module's test data directory fixture.
        json_data: JSON test file parsing fixture.
        ncbi_file: NCBI dataset genome statistics JSON file.
        core_file: Core database genome statistics JSON file.
        output_file: Expected output when comparing both sources.

    """
    ncbi_path = data_dir / ncbi_file
    core_path = data_dir / core_file
    result = compare.compare_stats_files(ncbi_path, core_path)
    output_data = json_data(output_file)
    assert not DeepDiff(result, output_data)
