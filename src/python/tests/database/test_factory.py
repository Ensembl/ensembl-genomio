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
"""Unit testing of `ensembl.io.genomio.database.factory` module."""
# pylint: disable=too-many-positional-arguments

from pathlib import Path
from unittest.mock import call, Mock, patch

from deepdiff import DeepDiff
import pytest
from pytest import param
from _pytest.capture import CaptureFixture
from sqlalchemy.engine import make_url, URL

from ensembl.io.genomio.database import factory


_META = {
    "species.production_name": "dog",
    "species.division": "metazoa",
    "assembly.accession": "GCA_000111222.3",
    "BRC4.organism_abbrev": "brc_dog",
    "BRC4.component": "brc_db",
}


@patch("ensembl.io.genomio.database.factory.DBConnectionLite")
@pytest.mark.parametrize(
    "server_url, dbs, brc_mode, skip_keys, output",
    [
        param(URL.create("mysql"), [], False, False, [], id="No databases selected"),
        param(
            URL.create("mysql"),
            ["db1"],
            False,
            False,
            [
                {
                    "server": {"host": None, "user": None, "port": None, "password": None, "database": "db1"},
                    "production_name": "dog",
                    "species": "dog",
                    "division": "metazoa",
                    "accession": "GCA_000111222.3",
                    "release": "110",
                }
            ],
            id="Ensembl core database",
        ),
        param(
            URL.create("mysql"),
            ["db1", "db2"],
            True,
            False,
            [
                {
                    "server": {"host": None, "user": None, "port": None, "password": None, "database": "db1"},
                    "production_name": "dog",
                    "species": "brc_dog",
                    "division": "brc_db",
                    "accession": "GCA_000111222.3",
                    "release": "110",
                },
                {
                    "server": {"host": None, "user": None, "port": None, "password": None, "database": "db2"},
                    "production_name": "dog",
                    "species": "brc_dog",
                    "division": "brc_db",
                    "accession": "GCA_000111222.3",
                    "release": "110",
                },
            ],
            id="VEuPathDB core databases",
        ),
        param(
            URL.create("mysql"),
            ["db1", "db2"],
            True,
            True,
            [
                {
                    "server": {"host": None, "user": None, "port": None, "password": None, "database": "db1"},
                    "production_name": "dog",
                    "species": "dog",
                    "division": "all",
                    "accession": "GCA_000111222.3",
                    "release": "110",
                },
                {
                    "server": {"host": None, "user": None, "port": None, "password": None, "database": "db2"},
                    "production_name": "dog",
                    "species": "dog",
                    "division": "all",
                    "accession": "GCA_000111222.3",
                    "release": "110",
                },
            ],
            id="VEuPathDB core databases",
        ),
    ],
)
def test_format_db_data(
    mock_dbconn: Mock, server_url: URL, dbs: list[str], brc_mode: bool, skip_keys: bool, output: list[dict]
) -> None:
    """Tests the `factory.format_db_data()` function.

    Args:
        mock_dbconn: A mock of `ensembl.io.genomio.database.factory.DBConnectionLite` class.
        server_url: Server URL where all the databases are hosted.
        dbs: List of database names.
        brc_mode: BRC mode?
        skip_keys: Return `None` instead of the assigned value for "BRC4.*" meta keys.
        output: Expected list of dictionaries with metadata per database.
    """

    def _get_meta_value(meta_key: str) -> str | None:
        """Return empty string if "species.division" is requested in BRC mode, "Metazoa" otherwise."""
        if (meta_key == "species.division") and brc_mode:
            return ""
        if meta_key.startswith("BRC4.") and skip_keys:
            return None
        return _META[meta_key]

    dbconnection = Mock()
    dbconnection.get_meta_value.side_effect = _get_meta_value
    dbconnection.get_project_release.return_value = "110"
    mock_dbconn.return_value = dbconnection

    result = factory.format_db_data(server_url, dbs, brc_mode)
    assert not DeepDiff(result, output)
    if dbs:
        calls = [call("species.production_name"), call("species.division"), call("assembly.accession")]
        if brc_mode:
            calls += [call("BRC4.organism_abbrev"), call("BRC4.component")]
        dbconnection.get_meta_value.assert_has_calls(calls)
        dbconnection.get_project_release.assert_called()


@patch("ensembl.io.genomio.database.factory.format_db_data")
@patch("ensembl.io.genomio.database.factory.CoreServer")
@pytest.mark.parametrize(
    "use_db_file, output",
    [
        param(
            False,
            [{"database": "db1", "species": "dog"}, {"database": "db2", "species": "dog"}],
            id="Get metadata for all databases",
        ),
        param(True, [{"database": "db1", "species": "dog"}], id="Use file to filter databases"),
    ],
)
def test_get_core_dbs_metadata(
    mock_core_server: Mock,
    mock_format_db_data: Mock,
    data_dir: Path,
    use_db_file: bool,
    output: list[dict],
) -> None:
    """Tests the `factory.get_core_dbs_metadata()` function.

    Args:
        mock_core_server: A mock of `ensembl.io.genomio.database.factory.CoreServer` class.
        mock_format_db_data: A mock of `ensembl.io.genomio.database.factory.format_db_data` function.
        data_dir: Module's test data directory fixture.
        use_db_file: Use database file to filter databases.
        output: Expected list of dictionaries with some metadata for each selected database.
    """

    def _format_db_data(server_url: URL, dbs: list[str], brc_mode: bool = False) -> list[dict]:
        """Returns metadata from a list of databases."""
        _ = (server_url, brc_mode)  # Unused by mock
        return [{"database": db, "species": "dog"} for db in dbs]

    if use_db_file:
        mock_core_server.get_cores.return_value = ["db1"]
    else:
        mock_core_server.get_cores.return_value = ["db1", "db2"]
    mock_core_server.return_value = mock_core_server
    mock_format_db_data.side_effect = _format_db_data

    db_list = data_dir / "databases.txt" if use_db_file else None
    result = factory.get_core_dbs_metadata(URL.create("mysql"), db_list=db_list)
    assert not DeepDiff(result, output)


@pytest.mark.parametrize(
    "arg_list, expected",
    [
        param(
            ["--host", "localhost", "--port", "42", "--user", "me"],
            {
                "host": "localhost",
                "port": 42,
                "user": "me",
                "password": None,
                "url": make_url("mysql://me@localhost:42"),
                "prefix": "",
                "build": None,
                "release": None,
                "db_regex": "",
                "db_list": None,
                "brc_mode": False,
                "log_level": "WARNING",
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
                "--password",
                "secret",
                "--prefix",
                "pre",
                "--build",
                "70",
                "--release",
                "114",
                "--db_regex",
                "%_core_%",
                "--db_list",
                __file__,
                "--brc_mode",
            ],
            {
                "host": "localhost",
                "port": 42,
                "user": "me",
                "password": "secret",
                "url": make_url("mysql://me:secret@localhost:42"),
                "prefix": "pre",
                "build": 70,
                "release": 114,
                "db_regex": "%_core_%",
                "db_list": __file__,
                "brc_mode": True,
                "log_level": "WARNING",
            },
            id="New arg values",
        ),
    ],
)
def test_parse_args(arg_list: list[str], expected: dict) -> None:
    """Tests the `factory.parse_args()` function."""
    args = factory.parse_args(arg_list)
    if args.db_list:
        # DeepDiff is not able to compare two objects of Path type, so convert it to string
        setattr(args, "db_list", str(args.db_list))
    assert not DeepDiff(vars(args), expected)


@pytest.mark.parametrize(
    "arg_list, server_url, stdout",
    [
        (
            ["--host", "localhost", "--port", "42", "--user", "me"],
            make_url("mysql://me@localhost:42"),
            '{\n    "test": "output"\n}\n',
        ),
    ],
)
@patch("ensembl.io.genomio.database.factory.get_core_dbs_metadata")
def test_main(
    mock_get_core_dbs_metadata: Mock,
    capsys: CaptureFixture[str],
    arg_list: list[str],
    server_url: URL,
    stdout: str,
) -> None:
    """Tests the `factory.main()` function (entry point).

    Fixtures: capsys
    """
    mock_get_core_dbs_metadata.return_value = {"test": "output"}
    factory.main(arg_list)
    # Check that we have called the mocked function once with the expected parameters
    mock_get_core_dbs_metadata.assert_called_once_with(
        server_url=server_url, prefix="", build=None, version=None, db_regex="", db_list=None, brc_mode=False
    )
    # Check that the stdout is as expected
    captured = capsys.readouterr()
    assert captured.out == stdout
