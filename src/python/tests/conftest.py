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

from pathlib import Path

import pytest
from pytest import Config, FixtureRequest


TEST_DATA_DIR_PATH = Path(__file__).parent
FILES_DIR_PATH = TEST_DATA_DIR_PATH / "data" / "flatfiles"


@pytest.fixture(scope="package")
def dbs_dir() -> Path:
    """Returns the folder that contains the database test files."""
    return TEST_DATA_DIR_PATH / "databases"


@pytest.fixture(scope="package")
def files_dir() -> Path:
    """Returns the folder that contains the flat test files."""
    return FILES_DIR_PATH


@pytest.fixture(scope="package")
def manifest_dir() -> Path:
    """Returns the folder that contains the manifest data test files."""
    return FILES_DIR_PATH / "manifest_data"


@pytest.fixture(scope="module")
def datadir(request: FixtureRequest) -> Path:
    """Returns the path to the test data folder matching the test's name.

    Args:
        request: Fixture providing information of the requesting test function.

    """
    return Path(request.module.__file__).with_suffix("")


@pytest.fixture(scope="package")
def shared_datadir(pytestconfig: Config) -> Path:
    """Returns the path to the shared test data folder.

    Args:
        pytestconfig: Session-scoped fixture that returns the session's `pytest.Config` object.

    """
    return Path(pytestconfig.rootpath) / "src" / "python" / "tests" / "data"