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
"""Unit tests for `ensembl.io.genomio.database.analysis`."""

from unittest.mock import MagicMock, Mock, patch

import pytest
from sqlalchemy.engine import URL

from ensembl.io.genomio.database import analysis


@patch("ensembl.io.genomio.database.analysis.DBConnection")
def test_get_production_analysis_description_returns_analysis_data(mock_dbconnection: Mock) -> None:
    """Return the requested analysis description and associated web data."""
    result = {
        "description": "Repeat finding",
        "display_label": "Repeats",
        "displayable": 1,
        "web_data": "{\"url\": \"https://example.org\"}",
    }
    session = Mock()
    session.execute.return_value.mappings.return_value.one_or_none.return_value = result
    database = Mock()
    database.session_scope.return_value = MagicMock()
    database.session_scope.return_value.__enter__.return_value = session
    mock_dbconnection.return_value = database
    prod_url = URL.create("mysql", database="production")

    actual = analysis.get_production_analysis_description("repeat_finder", prod_url)

    assert actual == result
    mock_dbconnection.assert_called_once_with(prod_url, reflect=False)
    session.execute.assert_called_once()
    query, parameters = session.execute.call_args.args
    assert "LEFT JOIN web_data" in str(query)
    assert "WHERE logic_name = :logic_name" in str(query)
    assert parameters == {"logic_name": "repeat_finder"}


@patch("ensembl.io.genomio.database.analysis.DBConnection")
def test_get_production_analysis_description_preserves_null_columns(mock_dbconnection: Mock) -> None:
    """Keep nullable description and web data values as `None`."""
    result = {
        "description": None,
        "display_label": "Repeats",
        "displayable": 0,
        "web_data": None,
    }
    session = Mock()
    session.execute.return_value.mappings.return_value.one_or_none.return_value = result
    database = Mock()
    database.session_scope.return_value = MagicMock()
    database.session_scope.return_value.__enter__.return_value = session
    mock_dbconnection.return_value = database

    actual = analysis.get_production_analysis_description("repeat_finder", URL.create("mysql"))

    assert actual == result


@patch("ensembl.io.genomio.database.analysis.DBConnection")
def test_get_production_analysis_description_raises_when_analysis_is_missing(mock_dbconnection: Mock) -> None:
    """Raise a clear error when the logic name is absent from production."""
    session = Mock()
    session.execute.return_value.mappings.return_value.one_or_none.return_value = None
    database = Mock()
    database.session_scope.return_value = MagicMock()
    database.session_scope.return_value.__enter__.return_value = session
    mock_dbconnection.return_value = database

    with pytest.raises(ValueError, match="Analysis not found in production database: missing_analysis"):
        analysis.get_production_analysis_description("missing_analysis", URL.create("mysql"))
