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
"""Unit testing of `ensembl.io.genomio.schemas.json` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_json.py

"""

from contextlib import nullcontext as does_not_raise
from os import PathLike
from pathlib import Path
from typing import ContextManager, List

from jsonschema.exceptions import ValidationError
import pytest
from pytest import raises

from ensembl.io.genomio.schemas import json


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
def test_schema_factory(
    tmp_path: Path, data_dir: Path, metadata_types: List[str], output: List[PathLike]
) -> None:
    """Tests the `schema_factory()` method.

    Args:
        tmp_path: Test's unique temporary directory fixture.
        data_dir: Module's test data directory fixture.
        manifest_dir: Path to the folder with the manifest JSON file to check.
        metadata_types: Metadata types to extract from `manifest` as JSON files.
        output: Expected created files.

    """
    json.schema_factory(data_dir, metadata_types, tmp_path)
    for file_name in output:
        print(f"Check {file_name} in {tmp_path}")
        assert (tmp_path / file_name).exists()


@pytest.mark.parametrize(
    "json_file, json_schema, expected",
    [
        ("seq_region.json", "seq_region", does_not_raise()),
        ("functional_annotation.json", "functional_annotation_schema.json", does_not_raise()),
        ("seq_region.json", "functional_annotation", raises(ValidationError)),
    ],
)
def test_schema_validator(data_dir: Path, json_file: str, json_schema: str, expected: ContextManager) -> None:
    """Tests the `schema_validator()` method.

    Args:
        data_dir: Module's test data directory fixture.
        json_file: Path to the JSON file to check.
        json_schema: JSON schema to validate `json_file` against.
        expected: Context manager for the expected exception, i.e. the test will only pass if that
            exception is raised. Use `contextlib.nullcontext` if no exception is expected.

    """
    json_path = data_dir / json_file
    if Path(json_schema).suffix == ".json":
        json_schema = str(data_dir / json_schema)
    with expected:
        json.schema_validator(json_path, json_schema)
