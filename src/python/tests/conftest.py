# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Local directory-specific hook implementations.

Since this file is located at the root of all ensembl-genomio Python tests, every test in every subfolder
will have access to the plugins, hooks and fixtures defined here.

"""

import json
from pathlib import Path
from typing import Any, Callable

import pytest
from pytest import Config

from ensembl.io.genomio.utils import get_json


@pytest.fixture(scope="package")
def shared_data_dir(pytestconfig: Config) -> Path:
    """Returns the path to the shared test data folder.

    Args:
        pytestconfig: Session-scoped fixture that returns the session's `pytest.Config` object.

    """
    return pytestconfig.rootpath / "src/python/tests/data"


@pytest.fixture(name="json_data")
def fixture_json_data(data_dir: Path) -> Callable[[str], Any]:
    """Returns a JSON test object factory.

    Args:
        data_dir: Module's test data directory fixture.

    """

    def _json_data(file_name: str) -> Any:
        """Returns the parsed JSON object from the given JSON test file."""
        return get_json(data_dir / file_name)

    return _json_data


class MockResponse:
    """Mock a `requests` response."""

    def __init__(self, json_str: str) -> None:
        """Store the JSON text response.

        Args:
            json_str: Expected JSON test response.
        """
        self.text = json_str

    @staticmethod
    def raise_for_status() -> None:
        """Mock, never raise Exception here."""

    def json(self) -> dict:
        """Returns the data decoded from the JSON text."""
        return json.loads(self.text)


@pytest.fixture(name="mock_response")
def mock_response() -> Callable:
    """Fixture to mock a `requests` response."""

    def _mock_response(json_str: str) -> MockResponse:
        return MockResponse(json_str)

    return _mock_response
