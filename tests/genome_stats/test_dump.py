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
"""Unit testing of `ensembl.io.genomio.genome_stats.dump` module.

Typical usage example::
    $ pytest test_dump.py

"""

from dataclasses import dataclass
from pathlib import Path
from string import Template
from typing import Any, Callable, Dict, List, Tuple
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy.orm import Session
from sqlalchemy.sql.expression import ClauseElement

from ensembl.core.models import Gene, Transcript
from ensembl.io.genomio.genome_stats import dump
from ensembl.io.genomio.utils import get_json


@dataclass
class MockResult:
    """Mocker of `sqlalchemy.engine.Result` class."""

    rows: List

    def __iter__(self) -> Any:
        """Iterates over the elements in `rows` attribute."""
        yield from self.rows

    def one(self) -> Tuple:
        """Returns the first element in `rows` attribute."""
        return self.rows[0]


# Define templates for all expected queries
ATTRIB_COUNTS_QUERY = Template(
    "SELECT seq_region_attrib.value, count(*) AS count_1 "
    "FROM seq_region_attrib JOIN attrib_type ON attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id"
    " WHERE attrib_type.code = '${code}' GROUP BY seq_region_attrib.value"
)
BIOTYPES_QUERY = Template(
    "SELECT ${table}.biotype, count(*) AS count_1 FROM ${table} GROUP BY ${table}.biotype"
)
FEATURE_STATS_TOTAL_QUERY = Template("SELECT count(*) AS count_1 FROM ${table}")
FEATURE_STATS_NULL_QUERY = Template(
    "SELECT count(*) AS count_1 FROM ${table} WHERE ${table}.description IS NULL"
)
FEATURE_STATS_SOURCE_QUERY = Template(
    "SELECT count(*) AS count_1 FROM ${table} WHERE ${table}.description LIKE '%[Source:%'"
)


class MockSession(Session):
    """Mocker of `sqlalchemy.orm.Session` class that replaces its `execute()` method for testing."""

    # pylint: disable-next=too-many-return-statements
    def execute(self, statement: ClauseElement) -> MockResult:  # type: ignore[override]
        """Returns a `MockResult` object representing results of the statement execution.

        Args:
            statement: An executable statement.

        """
        compiled_statement = str(statement.compile(compile_kwargs={"literal_binds": True}))
        query = compiled_statement.replace("\n", "")
        # Expected statements for test_get_attrib_counts()
        if query == ATTRIB_COUNTS_QUERY.substitute(code="coord_system_tag"):
            return MockResult([["chromosome", 7], ["scaffold", 2]])
        if query in (
            ATTRIB_COUNTS_QUERY.substitute(code="sequence_location"),
            ATTRIB_COUNTS_QUERY.substitute(code="codon_table"),
        ):
            return MockResult([])
        # Expected statements for test_get_biotypes()
        if query == BIOTYPES_QUERY.substitute(table="gene"):
            return MockResult([["protein_coding", 40], ["pseudogene", 1]])
        if query == BIOTYPES_QUERY.substitute(table="transcript"):
            return MockResult([])
        # Expected statements for test_get_feature_stats()
        if query == FEATURE_STATS_TOTAL_QUERY.substitute(table="gene"):
            return MockResult([[10]])
        if query == FEATURE_STATS_NULL_QUERY.substitute(table="gene"):
            return MockResult([[2]])
        if query == FEATURE_STATS_SOURCE_QUERY.substitute(table="gene"):
            return MockResult([[1]])
        if query in (
            FEATURE_STATS_TOTAL_QUERY.substitute(table="transcript"),
            FEATURE_STATS_NULL_QUERY.substitute(table="transcript"),
            FEATURE_STATS_SOURCE_QUERY.substitute(table="transcript"),
        ):
            return MockResult([[0]])
        # Fail test if the provided statement is not expected
        raise RuntimeError(f"Could not mock results from statement: {statement}")


class TestStatsGenerator:
    """Tests for the `StatsGenerator` class."""

    stats_gen: dump.StatsGenerator
    genome_stats: Dict[str, Any]

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, data_dir: Path) -> None:
        """Loads the required fixtures and values as class attributes.

        Args:
            data_dir: Module's test data directory fixture.

        """
        type(self).stats_gen = dump.StatsGenerator(MockSession())
        type(self).genome_stats = get_json(data_dir / "genome_stats.json")

    @pytest.mark.parametrize(
        "stats, output",
        [
            ({}, {}),
            (
                {"coord_system": {"supercontig": 2, "scaffold": 3}},
                {"coord_system": {"supercontig": 2, "scaffold": 3}},
            ),
            ({"coord_system": {"supercontig": 5}}, {"coord_system": {"scaffold": 5}}),
        ],
    )
    @pytest.mark.dependency(name="fix_scaffolds")
    def test_fix_scaffolds(self, stats: Dict, output: Dict) -> None:
        """Tests the `StatsGenerator._fix_scaffolds()` static method.

        Args:
            stats: Input statistic dictionary.
            output: Expected output.

        """
        dump.StatsGenerator._fix_scaffolds(stats)  # pylint: disable=protected-access
        assert stats == output

    @pytest.mark.parametrize(
        "code, attribute",
        [
            ("coord_system_tag", "coord_system"),
            ("sequence_location", "locations"),
        ],
    )
    @pytest.mark.dependency(name="get_attrib_counts")
    def test_get_attrib_counts(self, code: str, attribute: str) -> None:
        """Tests the `StatsGenerator.get_attrib_counts()` method.

        Args:
            code: Core database attribute type code.
            attribute: Core database attribute name.

        """
        attrib_counts = self.stats_gen.get_attrib_counts(code)
        assert attrib_counts == self.genome_stats["assembly_stats"][attribute]

    @pytest.mark.parametrize(
        "table, table_name",
        [
            (Gene, "genes"),
            (Transcript, "transcripts"),
        ],
    )
    @pytest.mark.dependency(name="get_biotypes")
    def test_get_biotypes(self, table: Any, table_name: str) -> None:
        """Tests the `StatsGenerator.get_biotypes()` method.

        Args:
            table: Core database table model class.
            table_name: Core database table name.

        """
        biotypes = self.stats_gen.get_biotypes(table)
        assert biotypes == self.genome_stats["annotation_stats"][table_name]["biotypes"]

    @pytest.mark.parametrize(
        "table, table_name",
        [
            (Gene, "genes"),
            (Transcript, "transcripts"),
        ],
    )
    @pytest.mark.dependency(name="get_feature_stats", depends=["get_biotypes"])
    def test_get_feature_stats(self, table: Any, table_name: str) -> None:
        """Tests the `StatsGenerator.get_feature_stats()` method.

        Args:
            table: Core database table model class.
            table_name: Core database table name.

        """
        biotypes = self.stats_gen.get_feature_stats(table)
        assert biotypes == self.genome_stats["annotation_stats"][table_name]

    @pytest.mark.dependency(name="get_assembly_stats", depends=["fix_scaffolds", "get_attrib_counts"])
    def test_get_assembly_stats(self) -> None:
        """Tests the `StatsGenerator.get_assembly_stats()` method."""
        assembly_stats = self.stats_gen.get_assembly_stats()
        assert assembly_stats == self.genome_stats["assembly_stats"]

    @pytest.mark.dependency(name="get_annotation_stats", depends=["get_feature_stats"])
    def test_get_annotation_stats(self) -> None:
        """Tests the `StatsGenerator.get_annotation_stats()` method."""
        annotation_stats = self.stats_gen.get_annotation_stats()
        assert annotation_stats == self.genome_stats["annotation_stats"]

    @pytest.mark.dependency(name="get_genome_stats", depends=["get_assembly_stats", "get_annotation_stats"])
    def test_get_genome_stats(self) -> None:
        """Tests the `StatsGenerator.get_genome_stats()` method."""
        genome_stats = self.stats_gen.get_genome_stats()
        assert genome_stats == self.genome_stats


@pytest.mark.dependency(depends=["get_genome_stats"])
@patch("ensembl.io.genomio.genome_stats.dump.DBConnectionLite")
def test_dump_genome_stats(mock_dbconnection: MagicMock, json_data: Callable) -> None:
    """Tests the `dump_genome_stats()` method.

    Args:
        mock_dbconnection: Mock of DBConnectionLite to avoid needing an actual core database.
        json_data: JSON test file parsing fixture.

    """
    mock_dbconnection.return_value.session_scope.return_value = MockSession()
    stats = dump.dump_genome_stats("")
    expected_output = json_data("genome_stats.json")
    assert stats == expected_output
