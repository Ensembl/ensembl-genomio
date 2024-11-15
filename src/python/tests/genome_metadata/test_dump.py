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
"""Unit testing of `ensembl.io.genomio.genome_metadata.dump` module.

Typical usage example::
    $ pytest test_dump.py

"""
# pylint: disable=too-many-positional-arguments

from collections import namedtuple
from contextlib import nullcontext as does_not_raise
from typing import Any, ContextManager
from unittest.mock import Mock, patch

from deepdiff import DeepDiff
import pytest
from pytest import param
from _pytest.capture import CaptureFixture
from sqlalchemy.engine import make_url, URL

from ensembl.io.genomio.genome_metadata import dump
from ensembl.utils import StrPath


MetaRow = namedtuple("MetaRow", "meta_key meta_value")


@pytest.mark.parametrize(
    ("genome_metadata", "output", "expectation"),
    [
        pytest.param({"assembly": {"version": "1"}}, 1, does_not_raise(), id="Version is '1'"),
        pytest.param(
            {"assembly": {"accession": "GCA_00000001.1", "version": "a"}},
            1,
            does_not_raise(),
            id="Version is 'a', accession's version is 1",
        ),
        pytest.param(
            {"assembly": {"accession": "GCA_00000001.1"}},
            1,
            does_not_raise(),
            id="No version, accession's version is 1",
        ),
        pytest.param(
            {"assembly": {"accession": "GCA_00000001"}},
            0,
            pytest.raises(ValueError),
            id="No version, accession without version",
        ),
    ],
)
def test_check_assembly_version(
    genome_metadata: dict[str, Any], output: int, expectation: ContextManager
) -> None:
    """Tests the `dump.check_assembly_version()` method.

    Args:
        genome_metadata: Nested genome metadata key values.
        output: Expected assembly version.
        expectation: Context manager for the expected exception (if any).
    """
    with expectation:
        dump.check_assembly_version(genome_metadata)
        assert genome_metadata["assembly"]["version"] == output


@pytest.mark.parametrize(
    ("genome_metadata", "output", "expectation"),
    [
        pytest.param({}, {}, does_not_raise(), id="No 'genebuild' entry"),
        pytest.param(
            {"genebuild": {"version": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="Version is 'v1', no ID",
        ),
        pytest.param(
            {"genebuild": {"version": "v1", "id": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="Version is 'v1', ID dropped",
        ),
        pytest.param(
            {"genebuild": {"id": "v1"}},
            {"genebuild": {"version": "v1"}},
            does_not_raise(),
            id="No version, ID moved to version",
        ),
        pytest.param({"genebuild": {}}, {}, pytest.raises(ValueError), id="No version or ID"),
    ],
)
def test_check_genebuild_version(
    genome_metadata: dict[str, Any], output: dict[str, Any], expectation: ContextManager
) -> None:
    """Tests the `dump.check_genebuild_version()` method.

    Args:
        genome_metadata: Nested genome metadata key values.
        output: Expected change in the genome metadata dictionary.
        expectation: Context manager for the expected exception (if any).
    """
    with expectation:
        dump.check_genebuild_version(genome_metadata)
        assert not DeepDiff(genome_metadata, output)


@patch("ensembl.io.genomio.genome_metadata.dump.check_genebuild_version", Mock())
@patch("ensembl.io.genomio.genome_metadata.dump.check_assembly_version", Mock())
@pytest.mark.parametrize(
    ("genome_metadata", "output", "metafilter", "meta_update"),
    [
        pytest.param(
            {"species": {"taxonomy_id": "5485"}},
            {"species": {"taxonomy_id": 5485}},
            None,
            False,
            id="Meta matches, no filter, allow meta update",
        ),
        pytest.param(
            {"species": {"taxonomy_id": "5485"}},
            {"species": {"taxonomy_id": 5485}},
            None,
            True,
            id="Meta matches, no filter, perform meta update",
        ),
        pytest.param(
            {"genebuild": {"new_key": "_"}}, {"genebuild": {}}, None, False, id="Filters on '_' value"
        ),
        pytest.param({"BRC5": "new_value"}, {}, None, False, id="BRC5 new value"),
        pytest.param(
            {"meta": "key", "species": {"alias": "woof"}},
            {"species": {"alias": "woof"}},
            None,
            False,
            id="Test alias",
        ),
        pytest.param(
            {"added_seq": {"region_name": [1, 2]}},
            {"added_seq": {"region_name": ["1", "2"]}},
            None,
            False,
            id="Added seq region_name",
        ),
        pytest.param({}, {}, None, False, id="BRC5 new value"),
        pytest.param(
            {
                "species": {
                    "display_name": "Honeybee",
                    "annotation_source": "Ensembl",
                    "production_name": "apis_melifera_gca123v1",
                    "scientific_name": "apis_melifera",
                    "taxonomy_id": "70921",
                }
            },
            {
                "species": {
                    "display_name": "Honeybee",
                    "production_name": "apis_melifera_gca123v1",
                    "taxonomy_id": "70921",
                }
            },
            {"species": {"display_name": "str", "production_name": "str", "taxonomy_id": "int"}},
            False,
            id="Filter via input meta JSON",
        ),
        pytest.param(
            {
                "annotation": {
                    "provider_name": "ENA",
                },
                "assembly": {"accession": "GCA_000111222.3", "version": "1"},
                "genebuild": {"method": "import", "version": "1"},
            },
            {"assembly": {"version": "1"}, "genebuild": {"method": "import", "version": "1"}},
            {"assembly": {"version": "str"}, "genebuild": {"version": "str", "method": "str"}},
            False,
            id="Asm + Genebuild version filter",
        ),
        pytest.param(
            {
                "annotation": {
                    "provider_name": "ENA",
                },
                "assembly": {"accession": "GCA_000111222.3", "version": "1"},
                "genebuild": {"method": "import", "version": "1"},
            },
            {"genebuild": {"method": "import"}},
            {"genebuild": {"method": "str"}},
            True,
            id="Only genebuild method, perform meta update",
        ),
    ],
)
def test_filter_genome_meta(
    genome_metadata: dict[str, Any],
    output: dict[str, Any],
    metafilter: StrPath,
    meta_update: bool,
) -> None:
    """Tests the `dump.filter_genome_meta()` method.

    Args:
        genome_metadata: Nested genome metadata key values.
        output: Expected change in the genome metadata dictionary.
        metafilter: Type evaluated meta filter.
        meta_update: Permit meta updating.
    """
    result = dump.filter_genome_meta(genome_metadata, metafilter, meta_update)
    assert not DeepDiff(result, output)


@pytest.mark.parametrize(
    ("meta_dict", "expected_dict"),
    [
        pytest.param(
            {"k1": {"sk1": "str"}, "k2": {"sk2": "float"}, "k3": {"sk3": "int"}},
            "{'k1': {'sk1': <class 'str'>}, 'k2': {'sk2': <class 'float'>}, 'k3': {'sk3': <class 'int'>}}",
            id="Filter conversion",
        ),
    ],
)
def test_convert_dict(meta_dict: dict, expected_dict: dict) -> None:
    """Tests the `dump.convert_dict()` method.

    Args:
        meta_dict: Dict containing string based meta 'subkey' value pairs.
        expected_dict: Dict with converted 'subkey' class types.
    """
    convert_dict = dump.convert_dict(meta_dict)
    string_convert = str(convert_dict)
    assert not DeepDiff(string_convert, expected_dict)


@patch("sqlalchemy.engine.Result")
@patch("sqlalchemy.orm.Session")
@pytest.mark.parametrize(
    ("db_name", "meta_data", "output", "expectation"),
    [
        pytest.param(None, [], {}, does_not_raise(), id="Empty meta table"),
        pytest.param(
            "test_dbname_core_110_1",
            [],
            {"database": {"name": "test_dbname_core_110_1"}},
            does_not_raise(),
            id="db_name append, empty meta table",
        ),
        pytest.param(
            None,
            [
                [MetaRow("sample", "gene1")],
                [MetaRow("species.name", "dog")],
                [MetaRow("species.synonym", "puppy")],
            ],
            {"sample": "gene1", "species": {"name": "dog", "synonym": "puppy"}},
            does_not_raise(),
            id="Meta table with simple values",
        ),
        pytest.param(
            None,
            [
                [MetaRow("sample", "gene1")],
                [MetaRow("sample", "gene2")],
                [MetaRow("species.synonym", "dog")],
                [MetaRow("species.synonym", "puppy")],
            ],
            {"sample": ["gene1", "gene2"], "species": {"synonym": ["dog", "puppy"]}},
            does_not_raise(),
            id="Meta table with lists",
        ),
        pytest.param(
            None,
            [[MetaRow("species", "dog")], [MetaRow("species.synonym", "puppy")]],
            {},
            pytest.raises(ValueError),
            id="'species' and 'species.synonym' meta keys",
        ),
        pytest.param(
            "test_dbname_core_110_1",
            [
                [MetaRow("assembly.accession", "GCA_000111222.3")],
                [MetaRow("species.annotation_source", "Community")],
                [MetaRow("species.production_name", "genus_species_gca000111222v3cm")],
            ],
            {
                "assembly": {"accession": "GCA_000111222.3"},
                "database": {"name": "test_dbname_core_110_1"},
                "species": {
                    "annotation_source": "Community",
                    "production_name": "genus_species_gca000111222v3cm",
                },
            },
            does_not_raise(),
            id="db_name append to meta",
        ),
    ],
)
def test_get_genome_metadata(
    mock_session: Mock,
    mock_result: Mock,
    db_name: str | None,
    meta_data: list[MetaRow],
    output: dict[str, Any],
    expectation: ContextManager,
) -> None:
    """Tests the `dump.get_genome_metadata()` method.

    Args:
        mock_session: A mock of `sqlalchemy.orm.Session()` class.
        db_name: Target core database name.
        meta_data: `meta` table content in a list of named tuples.
        output: Expected genome metadata dictionary.
        expectation: Context manager for the expected exception (if any).
    """
    mock_result.unique.return_value = mock_result
    mock_result.all.return_value = meta_data
    mock_session.execute.return_value = mock_result
    with expectation:
        result = dump.get_genome_metadata(mock_session, db_name)
        assert not DeepDiff(result, output)


@pytest.mark.parametrize(
    "arg_list, expected",
    [
        param(
            ["--host", "localhost", "--port", "42", "--user", "me", "--database", "test_db"],
            {
                "host": "localhost",
                "port": 42,
                "user": "me",
                "password": None,
                "url": make_url("mysql://me@localhost:42/test_db"),
                "database": "test_db",
                "metafilter": None,
                "meta_update": False,
                "append_db": False,
                "log_file": None,
                "log_level": "WARNING",
                "log_file_level": "DEBUG",
            },
            id="Default args",
        ),
        param(
            [
                "--host",
                "localhost",
                "--port",
                "42",
                "--user",
                "me",
                "--database",
                "test_db",
                "--metafilter",
                f"{__file__}",
                "--append_db",
            ],
            {
                "host": "localhost",
                "port": 42,
                "user": "me",
                "password": None,
                "url": make_url("mysql://me@localhost:42/test_db"),
                "database": "test_db",
                "metafilter": __file__,
                "meta_update": False,
                "append_db": True,
                "log_file": None,
                "log_level": "WARNING",
                "log_file_level": "DEBUG",
            },
            id="Filter, non-default args",
        ),
    ],
)
def test_parse_args(arg_list: list[str], expected: dict) -> None:
    """Tests the `dump.parse_args()` function."""
    # pylint: disable=too-many-positional-arguments
    args = dump.parse_args(arg_list)
    if args.metafilter:
        # DeepDiff is not able to compare two objects of Path type, so convert it to string
        setattr(args, "metafilter", str(args.metafilter))
    assert not DeepDiff(vars(args), expected)


@pytest.mark.parametrize(
    ("arg_list", "db_url", "metafilter", "meta_update", "append_db", "stdout"),
    [
        param(
            [
                "--host",
                "localhost",
                "--port",
                "42",
                "--user",
                "me",
                "--database",
                "test_dbname_core_110_1",
                "--append_db",
            ],
            make_url("mysql://me@localhost:42/test_dbname_core_110_1"),
            None,
            False,
            True,
            '{\n  "database": {\n    "name": "test_dbname_core_110_1"\n  }\n}\n',
            id="Call main and append_db",
        ),
    ],
)
@patch("ensembl.io.genomio.genome_metadata.dump.metadata_dump_setup")
def test_main(
    mock_metadata_dump_setup: Mock,
    capsys: CaptureFixture[str],
    arg_list: list[str],
    db_url: URL,
    metafilter: StrPath,
    meta_update: bool,
    append_db: bool,
    stdout: str,
) -> None:
    """Tests the `dump.main()` function (entry point).

    Fixtures: capsys
    """
    # pylint: disable=too-many-positional-arguments
    mock_metadata_dump_setup.return_value = {"database": {"name": "test_dbname_core_110_1"}}
    dump.main(arg_list)
    # Check that we have called the mocked function once with the expected parameters
    mock_metadata_dump_setup.assert_called_once_with(
        db_url=db_url, input_filter=metafilter, meta_update=meta_update, append_db=append_db
    )
    # Check that the stdout is as expected
    captured = capsys.readouterr()
    assert captured.out == stdout
