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

test_data_dir_path = Path(__file__).parent
root_dir_path = test_data_dir_path.parents[2]
files_dir_path = test_data_dir_path / "flatfiles"


@pytest.fixture(scope="package")
def dbs_dir():
    """Dir with database test files."""
    return test_data_dir_path / "databases"


@pytest.fixture(scope="package")
def files_dir():
    """Dir with flat test files."""
    return files_dir_path


@pytest.fixture(scope="package")
def manifest_dir():
    """Dir with manifest data test files."""
    return files_dir_path / "manifest_data"


@pytest.fixture(scope="package")
def schemas_dir():
    """Dir with schema test files."""
    return root_dir_path / "schemas"
