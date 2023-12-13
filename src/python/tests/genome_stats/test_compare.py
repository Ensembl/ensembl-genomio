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

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_compare.py

"""
# pylint: disable=eval-used,protected-access

from pathlib import Path
from typing import Any, Dict

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_stats import compare
from ensembl.io.genomio.utils import get_json


class TestCompare:
    """Tests `genome_stats.compare` methods."""

    ncbi_annot_stats: Dict[str, Any]
    ncbi_unannot_stats: Dict[str, Any]
    core_annot_stats: Dict[str, Any]
    core_unannot_stats: Dict[str, Any]
    output_empty: Dict[str, Any]
    output_unannot_diff: Dict[str, Any]
    output_annot_same: Dict[str, Any]
    output_annot_diff: Dict[str, Any]

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, genome_stats_dir: Path) -> None:
        """Loads necessary fixtures and values as class attributes.

        Args:
            genome_stats_dir: Path to the genome stats test files (fixture).

        """
        type(self).ncbi_annot_stats = get_json(genome_stats_dir / "ncbi_annotated.json")
        type(self).ncbi_unannot_stats = get_json(genome_stats_dir / "ncbi_unannotated.json")
        type(self).core_annot_stats = get_json(genome_stats_dir / "core_annotated.json")
        type(self).core_unannot_stats = get_json(genome_stats_dir / "core_unannotated.json")
        type(self).output_empty = get_json(genome_stats_dir / "output_empty.json")
        type(self).output_unannot_diff = get_json(genome_stats_dir / "output_unannotated.json")
        type(self).output_annot_same = get_json(genome_stats_dir / "output_annotated_same.json")
        type(self).output_annot_diff = get_json(genome_stats_dir / "output_annotated_diff.json")

    @pytest.mark.dependency(name="test_compare_dicts", scope="class")
    @pytest.mark.parametrize(
        "ncbi, core, output",
        [
            pytest.param({}, {}, {"same": {}, "different": {}}, id="empty_dicts"),
            pytest.param({"a": 0}, {"a": 0}, {"same": {}, "different": {}}, id="same_dicts_zero_values"),
            pytest.param(
                {"a": 3}, {"a": 3}, {"same": {"a": 3}, "different": {}}, id="same_dicts_non_zero_values"
            ),
            pytest.param(
                {"a": 3},
                {"a": 5},
                {"same": {}, "different": {"a": {"ncbi": 3, "core": 5, "diff": 2}}},
                id="different_dicts",
            ),
            pytest.param(
                {"a": 3, "b": 5},
                {"a": 3, "b": 3},
                {"same": {"a": 3}, "different": {"b": {"ncbi": 5, "core": 3, "diff": -2}}},
                id="partially_similar_dicts",
            ),
        ],
    )
    def test_compare_dicts(self, ncbi: Dict[str, int], core: Dict[str, int], output: Dict[str, Dict]) -> None:
        """Tests the `compare._compare_dicts()` method.

        Args:
            ncbi: NCBI dataset statistics in key-value pairs.
            core: Core database statistics in key-value pairs.
            output: Expected output when comparing both dictionaries.

        """
        result = compare._compare_dicts(ncbi, core)
        assert not DeepDiff(result, output)

    @pytest.mark.dependency(name="test_compare_assembly", depends=["test_compare_dicts"], scope="class")
    @pytest.mark.parametrize(
        "ncbi, core, output",
        [
            pytest.param("{}", "{}", "self.output_empty['assembly_diff']", id="empty_stats"),
            pytest.param(
                "self.ncbi_unannot_stats['reports'][0]",
                "self.core_unannot_stats['assembly_stats']",
                "self.output_unannot_diff['assembly_diff']",
                id="contig_assembly",
            ),
            pytest.param(
                "self.ncbi_annot_stats['reports'][0]",
                "self.core_annot_stats['assembly_stats']",
                "self.output_annot_same['assembly_diff']",
                id="chromosome_assembly",
            ),
        ],
    )
    def test_compare_assembly(self, ncbi: str, core: str, output: str) -> None:
        """Tests the `compare.compare_assembly()` method with simple data.

        Args:
            ncbi: NCBI dataset assembly statistics.
            core: Core database assembly statistics.
            output: Expected output when comparing both sources.

        """
        result = compare.compare_assembly(eval(ncbi), eval(core))
        assert not DeepDiff(result, eval(output))

    @pytest.mark.dependency(name="test_compare_annotation", depends=["test_compare_dicts"], scope="class")
    @pytest.mark.parametrize(
        "ncbi, core, output",
        [
            pytest.param("{}", "{}", "self.output_empty['assembly_diff']", id="empty_stats"),
            pytest.param(
                "self.ncbi_annot_stats['reports'][0]['annotation_info']['stats']['gene_counts']",
                "self.core_annot_stats['annotation_stats']",
                "self.output_annot_same['annotation_diff']",
                id="annotated_assembly",
            ),
        ],
    )
    def test_compare_annotation(self, ncbi: str, core: str, output: str) -> None:
        """Tests the `compare.compare_annotation()` method when the assembly has no annotation.

        Args:
            ncbi: NCBI dataset annotation statistics.
            core: Core database annotation statistics.
            output: Expected output when comparing both sources.

        """
        result = compare.compare_annotation(eval(ncbi), eval(core))
        assert not DeepDiff(result, eval(output))

    @pytest.mark.dependency(
        name="test_compare_stats", depends=["test_compare_assembly", "test_compare_annotation"], scope="class"
    )
    @pytest.mark.parametrize(
        "ncbi, core, output",
        [
            pytest.param("{}", "{}", "self.output_empty", id="empty_stats"),
            pytest.param(
                "self.ncbi_annot_stats['reports'][0]",
                "self.core_annot_stats",
                "self.output_annot_same",
                id="annotated_assembly",
            ),
        ],
    )
    def test_compare_stats(self, ncbi: str, core: str, output: str) -> None:
        """Tests the `compare.compare_stats()` method.

        Args:
            ncbi: NCBI dataset genome statistics.
            core: Core database genome statistics.
            output: Expected output when comparing both sources.

        """
        result = compare.compare_stats(eval(ncbi), eval(core))
        assert not DeepDiff(result, eval(output))

    @pytest.mark.dependency(name="test_compare_stats_files", depends=["test_compare_stats"], scope="class")
    @pytest.mark.parametrize(
        "ncbi_file, core_file, output",
        [
            pytest.param(
                "empty_file.txt", "core_unannotated.json", "self.output_empty", id="empty_ncbi_file"
            ),
            pytest.param(
                "ncbi_empty.json", "core_unannotated.json", "self.output_empty", id="wrong_ncbi_format"
            ),
            pytest.param(
                "ncbi_no_reports.json", "core_annotated.json", "self.output_annot_diff", id="no_ncbi_reports"
            ),
            pytest.param(
                "ncbi_annotated.json", "core_annotated.json", "self.output_annot_same", id="annotated_genome"
            ),
        ],
    )
    def test_compare_stats_files(
        self, genome_stats_dir: Path, ncbi_file: str, core_file: str, output: str
    ) -> None:
        """Tests the `compare.compare_stats_files()` method.

        Args:
            genome_stats_dir: Path to the genome stats test files (fixture).
            ncbi_file: NCBI dataset genome statistics JSON file.
            core_file: Core database genome statistics JSON file.
            output: Expected output when comparing both sources.

        """
        result = compare.compare_stats_files(genome_stats_dir / ncbi_file, genome_stats_dir / core_file)
        assert not DeepDiff(result, eval(output))
