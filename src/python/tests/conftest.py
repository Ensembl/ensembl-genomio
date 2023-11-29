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

TEST_DATA_DIR_PATH = Path(__file__).parent
ROOT_DIR_PATH = TEST_DATA_DIR_PATH.parents[2]
FILES_DIR_PATH = TEST_DATA_DIR_PATH / "flatfiles"


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


@pytest.fixture(scope="package")
def schemas_dir() -> Path:
    """Returns the folder that contains the manifest data test files."""
    return FILES_DIR_PATH / "schemas"
