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
"""Unit testing of `ensembl.io.genomio.assembly.download` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest ./assembly/test_download.py

"""

from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from contextlib import nullcontext as does_not_raise
from typing import ContextManager

from typing import Dict

import pytest
from ftplib import FTP


from ensembl.io.genomio.assembly.download import FTPConnection, FTPConnectionError
from ensembl.io.genomio.assembly.download import get_checksums, md5_files


# class TestDownloadAssembly:
#     """Tests class for `download_assembly` module."""


@pytest.mark.parametrize(
    "ftp_url, sub_dir, expectation",
    [
        ("ftp.ncbi.nlm.nih.gov", "genomes/all/GCA/017/607/445", does_not_raise()),  # normal case
        (
            "ftp.ncbi.nlm.nih.gov",
            "genomes/all/GCA/017/607/445",
            pytest.raises(FTPConnectionError),
        ),  # abnormal case, should run, but fail assert
        ("ftp.ncbi.nlm.nih.gov", None, pytest.raises(FTPConnectionError)),  # bad case no subdir path
        (
            "ftp.remote.fake.gov",
            "genomes/all/GCA/017/607/445",
            pytest.raises(FTPConnectionError),
        ),  # bad case incorrect ftp_url
        (
            "ftp.ncbi.nlm.nih.gov",
            "genomes/GCA/017/607/445",
            pytest.raises(FTPConnectionError),
        ),  # bad case malformed subdir path
    ],
)
@patch("ensembl.io.genomio.assembly.download.FTP", autospec=True)
def test_ftp_connection(
    mock_ftp_constructor: Mock,
    ftp_url: str,
    sub_dir: str,
    expectation: ContextManager,
):
    """Tests the FTPConnection class method 'establish_ftp()'
    Args:
        mock_ftp_constructor:
        ftp_url:
        sub_dir:
    """
    mock_ftp = mock_ftp_constructor.return_value
    mock_ftp.pwd = MagicMock(return_value="ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/607/445/")
    connected_ftp = FTPConnection.establish_ftp(mock_ftp, ftp_url, sub_dir)
    listing = connected_ftp.pwd()
    assert listing == "ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/607/445/"
    connected_ftp.pwd.assert_called_once()

#################
@pytest.mark.dependency(name="test_checksums")
@pytest.mark.parametrize(
    "checksum_file, checksum, expectation",
    [
        (Path("md5checksums.txt"), "40df91d5c40cb55621c4c92201da6834", does_not_raise()),  # Normal case
        (
            Path("malformed_md5_checksums.txt"),
            None,
            pytest.raises(ValueError),
        ),  # Bad md5sum file (tab separated, expects two \w)
        (Path("dir_name"), None, does_not_raise()),  # pass a dir instead of checksum
    ],
)
def test_checksums(data_dir: Path, checksum_file: Path, checksum: str, expectation: ContextManager) -> None:
    """Tests the 'download.get_checksums() function"""

    with expectation:
        md5_input_path = data_dir / checksum_file
        obtained_checksums = get_checksums(md5_input_path)
        assert obtained_checksums.get("annotation_hashes.txt", None) == checksum

#################
@pytest.mark.dependency(name="test_md5_files")
@pytest.mark.parametrize(
    "md5file, checksum_bool",
    [
        ("md5checksums.txt", True),  # Normal case
        ("wrong_md5_checksums.txt", False),  # Bad md5sum file (in correct md5)
        (None, True),  # No md5file specified, resort to default
        (Path("*"), False),  # Pass checksum_file as incompatible os path '*'
        ("missingfile_md5.txt", False),  # Pass md5_checksum file with reference to missing file
    ],
)
@pytest.mark.dependency(depends=["test_checksums"])
def test_md5_files(data_dir: Path, md5file: Path, checksum_bool: bool) -> None:
    """Tests the md5_files() function"""
    return_bool_on_md5files = md5_files(data_dir, md5file)
    assert return_bool_on_md5files == checksum_bool

#################


