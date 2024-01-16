#!/usr/bin/env python
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

from ensembl.io.genomio.schemas import json


class TestJSONSchemas:
    """Tests for the schemas modules."""

    test_data_dir: Path
    tmp_dir: Path

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, tmp_dir: Path, schemas_dir: Path):
        """Loads necessary fixtures and values as class attributes."""
        type(self).test_data_dir = schemas_dir
        type(self).tmp_dir = tmp_dir

    @pytest.mark.parametrize(
        "metadata_types, output",
        [
            (["new_metadata"], ["manifest.json"]),
            (
                ["functional_annotation", "seq_region"],
                ["manifest.json", "functional_annotation_category.json", "seq_region.json"],
            ),
        ],
    )
    def test_schema_factory(self, metadata_types: List[str], output: List[PathLike]) -> None:
        """Tests `ensembl.io.genomio.schemas.json.schema_factory()` method.

        Args:
            manifest_dir: Path to the folder with the manifest JSON file to check.
            metadata_types: Metadata types to extract from `manifest` as JSON files.
            output: Expected created files.

        """
        json.schema_factory(self.test_data_dir, metadata_types, self.tmp_dir)
        for file_name in output:
            print(f"Check {file_name} in {self.tmp_dir}")
            assert (self.tmp_dir / file_name).exists()

    @pytest.mark.parametrize(
        "json_file, json_schema, expected",
        [
            ("seq_region.json", "seq_region", does_not_raise()),
            ("functional_annotation.json", "functional_annotation_schema.json", does_not_raise()),
            ("seq_region.json", "functional_annotation", raises(ValidationError)),
        ],
    )
    def test_schema_validator(self, json_file: str, json_schema: str, expected: ContextManager) -> None:
        """Tests `ensembl.io.genomio.schemas.json.schema_validator()` method.

        Args:
            json_file: Path to the JSON file to check.
            json_schema: JSON schema to validate `json_file` against.
            expected: Context manager for the expected exception, i.e. the test will only pass if that
                exception is raised. Use :class:`~contextlib.nullcontext` if no exception is expected.

        """
        json_path = self.test_data_dir / json_file
        if Path(json_schema).suffix == ".json":
            json_schema = str(self.test_data_dir / json_schema)
        with expected:
            json.schema_validator(json_path, json_schema)