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
"""Unit testing of :mod:`ensembl.io.genomio.schemas` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_schemas.py

"""

from contextlib import nullcontext as does_not_raise
from os import PathLike
from pathlib import Path
from typing import ContextManager, List

from jsonschema.exceptions import ValidationError
import pytest
from pytest import raises

from ensembl.io.genomio import schemas


class TestSchemas:
    @pytest.fixture(scope="class", autouse=True)
    def setup(self, tmp_dir: Path):
        """Loads necessary fixtures and values as class attributes."""
        type(self).test_data_dir = pytest.files_dir / "schemas"
        type(self).tmp_dir = tmp_dir

    @pytest.mark.parametrize(
        "metadata_types, output",
        [
            (["new_metadata"], ["manifest.json"]),
            (
                ["functional_annotation", "seq_region"],
                ["manifest.json", "functional_annotation_test.json", "seq_region.json"],
            ),
        ],
    )
    def test_json_schema_factory(self, metadata_types: List[str], output: List[PathLike]) -> None:
        """Tests :meth:`schemas.json_schema_factory()` method.

        Args:
            manifest_dir: Path to the folder with the manifest JSON file to check.
            metadata_types: Metadata types to extract from `manifest` as JSON files.
            output: Expected created files.

        """
        schemas.json_schema_factory(self.test_data_dir, metadata_types, self.tmp_dir)
        for file_name in output:
            assert (self.tmp_dir / file_name).exists()

    @pytest.mark.parametrize(
        "json_file, json_schema, expected",
        [
            ("seq_region.json", "seq_region_schema.json", does_not_raise()),
            ("seq_region.json", "functional_annotation_schema.json", raises(ValidationError)),
        ],
    )
    def test_validate_json_schema(self, json_file: str, json_schema: str, expected: ContextManager) -> None:
        """Tests :meth:`schemas.validate_json_schema()` method.

        Args:
            json_file: Path to the JSON file to check.
            json_schema: JSON schema to validate `json_file` against.
            expected: Context manager for the expected exception, i.e. the test will only pass if that
                exception is raised. Use :class:`~contextlib.nullcontext` if no exception is expected.

        """
        with expected:
            schemas.validate_json_schema(self.test_data_dir / json_file, self.test_data_dir / json_schema)
