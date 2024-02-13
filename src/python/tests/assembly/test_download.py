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

import filecmp
import logging
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from contextlib import nullcontext as does_not_raise
from typing import ContextManager

from typing import Dict

from ftplib import FTP, error_reply as ftp_error_reply
import pytest


from ensembl.io.genomio.assembly.download import FTPConnection, FTPConnectionError
from ensembl.io.genomio.assembly.download import get_checksums, md5_files, _download_file


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
# @pytest.mark.dependency(name="test_download_file")
# @pytest.mark.parametrize(
#     "ftp_url, sub_dir, ftp_file, md5_sums",
#     [
#         (
#             "ftp.ncbi.nlm.nih.gov",
#             "genomes/all/GCA/017/607/445",
#             "GCA_017607445.1_ASM1760744v1_assembly_report.txt",
#             dict([('GCA_017607445.1_ASM1760744v1_assembly_report.txt','a03f39d1de753fcd380bf0775d5205d0')]),
#         ),
#         (
#             "ftp.ncbi.nlm.nih.gov",
#             "genomes/all/GCA/017/607/445",
#             "non-existing-file.txt",
#             dict([('GCA_017607445.1_ASM1760744v1_assembly_report.txt','a03f39d1de753fcd380bf0775d5205d0')]),
#         ),
#     ],
# )
    

# class TestDownloadFile(unittest.TestCase):

#     @classmethod
#     def setUpClass(cls):
#         # mock return value
#         cls.patcher = patch('ensembl.io.genomio.assembly.download.FTP', return_value='226 Transfer complete.', autospec=True)

#         # start the patch
#         cls.patcher.start()

#         # stop after all tests
#         cls.addClassCleanup(cls.patcher.stop)


@pytest.mark.parametrize(
    "ftp_file, md5_sums, expectation",
    [
        pytest.param(
            "test_ftp_file.txt",
            {"test_ftp_file.txt": "e98b980b442fdb2a21877dcc55e11848"},
            does_not_raise(),
            id="Download OK file",
        ),
        pytest.param(
            "non-existing-file.txt",
            {"test_ftp_file.txt": "e98b980b442fdb2a21877dcc55e11848"},
            pytest.raises(FileNotFoundError),
            id="Can't download non-existing file",
        ),
        pytest.param(
            "non-existing-file.txt",
            {"non-existing-file.txt": "e98b980b442fdb2a21877dcc55e11848"},
            pytest.raises(ftp_error_reply),
            id="File has md5sum but does not exist",
        ),
    ],
)
@patch("ftplib.FTP")
def test_download_single_file(
    mock_ftp: FTP, data_dir: Path, tmp_dir: Path, ftp_file: str, md5_sums: dict, expectation: ContextManager
) -> None:
    """Tests the private function _download_single_file."""

    data_file = data_dir / ftp_file
    retr_file = tmp_dir / ftp_file

    def mock_retr_binary(command: str, callback: object):
        logging.info(f"Faking the download of {command}")
        try:
            with data_file.open("rb") as data_fh:
                callback(data_fh.read())
        except OSError as err:
            raise ftp_error_reply from err

    mock_ftp.retrbinary = MagicMock(side_effect=mock_retr_binary)
    with expectation:
        _download_file(mock_ftp, ftp_file, md5_sums, tmp_dir)
        assert filecmp.cmp(data_file, retr_file)


# @pytest.mark.dependency(depends=["test_md5_files"])
# @patch("ensembl.io.genomio.assembly.download._download_file", autospec=True)
# def test_download_single_file(
#     mock_ftp_constructor: Mock, 
#     ftp_url: str, 
#     sub_dir: str,
#     ftp_file: str, 
#     md5_sums: Dict[str, str], 
#     data_dir: Path, 
# ) -> None:
#     """
#     """
#     mock_ftp = mock_ftp_constructor.return_value
#     mock_ftp.pwd = MagicMock(return_value="ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/607/445/")
#     connected_ftp = FTPConnection.establish_ftp(mock_ftp, ftp_url, sub_dir)

#     _download_file(connected_ftp, ftp_file, md5_sums, data_dir)
#     .is_file.assert_called_once_with(Path(data_dir, ftp_file))
#     connected_ftp.retrbinary.assert_called()
#     connected_ftp.retrbinary.assert_called()