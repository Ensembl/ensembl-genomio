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

# from typing import Dict

from ftplib import FTP, error_reply as ftp_error_reply
import pytest

from ensembl.io.genomio.assembly.download import establish_ftp
from ensembl.io.genomio.assembly.download import get_checksums, md5_files, _download_file
from ensembl.io.genomio.assembly.download import download_files, get_files_selection, retrieve_assembly_data
from ensembl.io.genomio.assembly.download import UnsupportedFormatError
from ensembl.io.genomio.assembly.download import FileDownloadError
from ensembl.io.genomio.assembly.download import FTPConnectionError


# class TestDownloadAssembly:
#     """Tests class for `download_assembly` module."""

#################
@pytest.mark.dependency(name="test_ftp_connection")
@pytest.mark.parametrize(
    "ftp_url, accession, expectation",
    [
        pytest.param(
            "ftp.ncbi.nlm.nih.gov", "GCA_017607445.1", does_not_raise(), id="Successful ftp connection"
        ),
        pytest.param(
            "", "GCA_017607445.1", pytest.raises(FTPConnectionError), id="Failed connection case bad url"
        ),
        pytest.param(
            "ftp.ncbi.nlm.nih.gov", "gcx_017607445.1", pytest.raises(UnsupportedFormatError), id="Failed connection improper INSDC accession"
        ),
    ],
)
@patch("ftplib.FTP")
def test_ftp_connection(
    mock_ftp: Mock,
    ftp_url: str,
    accession: str,
    expectation: ContextManager,
    ):
    """Tests the FTPConnection method 'establish_ftp()'
    Args:
        mock_ftp: Mock FTP Object 
        ftp_url: FTP url
        sub_dir: Sub directory path
        expectation: Context manager expected raise exception
    """

    def side_eff_conn(url:str):
        if not url:
            raise Exception()

    mock_ftp.connect.side_effect = side_eff_conn

    with expectation:
        connected_ftp = establish_ftp(mock_ftp, ftp_url, accession)
        connected_ftp.connect.assert_called_once_with(ftp_url)
        connected_ftp.login.assert_called_once()
        connected_ftp.cwd.assert_called_once_with("genomes/all/GCA/017/607/445")

#################
@pytest.mark.dependency(name="test_checksums")
@pytest.mark.parametrize(
    "checksum_file, checksum, expectation",
    [
        pytest.param(
            Path("md5checksums.txt"), "40df91d5c40cb55621c4c92201da6834", does_not_raise(), id="Normal case"
        ),
        pytest.param(
            Path("malformed_md5_checksums.txt"),
            None,
            pytest.raises(ValueError),
            id="Bad md5 file - tab separated",
        ),
        pytest.param(Path("dir_name"), None, does_not_raise(), id="Dir passed instead of file"),
    ],
)
def test_checksums(data_dir: Path, checksum_file: Path, checksum: str, expectation: ContextManager) -> None:
    """Tests the 'download.get_checksums() function
    
    Args:
        data_dir:
        checksum_file:
        checksum:
        expectation:
    """
    with expectation:
        md5_input_path = data_dir / checksum_file
        obtained_checksums = get_checksums(md5_input_path)
        assert obtained_checksums.get("annotation_hashes.txt", None) == checksum


#################
@pytest.mark.dependency(name="test_md5_files")
@pytest.mark.parametrize(
    "md5file, checksum_bool",
    [
        pytest.param("md5checksums.txt", True, id="Normal case"),
        pytest.param("wrong_md5_checksums.txt", False, id="Incorrect md5 checksum"),
        pytest.param(None, True, id="No md5file specified, resort to default"),
        pytest.param(Path("*"), False, id="Incompatible os path '*'"),
        pytest.param("missingfile_md5.txt", False, id="md5 checksum with ref of missing file"),
    ],
)
@pytest.mark.dependency(depends=["test_checksums"])
def test_md5_files(data_dir: Path, md5file: Path, checksum_bool: bool) -> None:
    """Tests the md5_files() function"""
    return_bool_on_md5files = md5_files(data_dir, md5file)
    assert return_bool_on_md5files == checksum_bool


#################
@pytest.mark.dependency(name="test_download_single_file")
@pytest.mark.parametrize(
    "ftp_file, md5_sums, max_redo, expectation",
    [
        pytest.param(
            "test_ftp_file.txt",
            {"test_ftp_file.txt": "e98b980b442fdb2a21877dcc55e11848"},
            2,
            does_not_raise(),
            id="Download OK file",
        ),
        pytest.param(
            "non-existing-file.txt",
            {"test_ftp_file.txt": "e98b980b442fdb2a21877dcc55e11848"},
            0,
            pytest.raises(FileNotFoundError),
            id="Can't download non-existing file",
        ),
        pytest.param(
            "non-existing-file.txt",
            {"non-existing-file.txt": "e98b980b442fdb2a21877dcc55e11848"},
            0,
            pytest.raises(ftp_error_reply),
            id="File has md5sum but does not exist",
        ),
        pytest.param(
            "test_ftp_file.txt",
            {"test_ftp_file.txt": ""},
            2,
            does_not_raise(),
            id="Max retry attempts non-zero",
        ),
        pytest.param(
            "test_ftp_file.txt",
            {"test_ftp_file.txt": "e98b980b442fdb2a21877XxxXxxXxxXx"},
            3,
            pytest.raises(FileDownloadError),
            id="Incorrect md5sum",
        ),
    ],
)
@patch("ftplib.FTP")
def test_download_single_file(
    mock_ftp: Mock, data_dir: Path, tmp_dir: Path, ftp_file: str, md5_sums: dict, max_redo: int, expectation: ContextManager
    ) -> None:
    """Tests the private function _download_file.
    
    Args:
        mock_ftp: Mock FTP object
        data_dir: Path to test data root dir
        tmp_dir: Temp dir created for test
        ftp_file: FTP file which to mock download
        md5_sums: FTP file and md5_sum value pair
        expectation: Context manager expected raise exception
    """

    data_file = data_dir / ftp_file
    retr_file = tmp_dir / ftp_file

    def mock_retr_binary(command: str, callback: object):
        logging.info(f"Faking the download of {command}")
        try:
            with data_file.open("rb") as data_fh:
                callback(data_fh.read())
        except OSError as err:
            raise ftp_error_reply from err

    mock_ftp.retrbinary.side_effect=mock_retr_binary

    with expectation:
        _download_file(mock_ftp, ftp_file, md5_sums, tmp_dir)
        assert filecmp.cmp(data_file, retr_file)

#################
@pytest.mark.dependency(name="test_download_all_files")
@pytest.mark.parametrize(
    "ftp_url, accession, md5, exception, max_redo",
    [
        pytest.param(
            "ftp.ncbi.nlm.nih.gov",
            "GCA_017607445.1",
            "md5checksums_ftp_source.txt",
            does_not_raise(),
            2,
            id="Normal case: Properly formatted GC[A/F] accession",
        ),
    ]
)
@pytest.mark.dependency(depends=["test_download_single_file"])
@patch("ftplib.FTP")
def test_download_all_files(mock_ftp: MagicMock, data_dir: Path, ftp_url: str,
                            md5: str, accession: str, exception: ContextManager, 
                            max_redo: int
                            ) -> None:
    
    data_file = data_dir / md5

    def side_eff_ftp_mlsd():
        
        mlsd_ret = []        
        files = data_dir.glob("*GCA_017607445.1*.gz")
        for file_path in files:
            ftp_file = str(file_path).split("/")[-1]
            mlsd_ret.append((str(ftp_file), ["fact_1","fact_2"]))

        return mlsd_ret

    mock_ftp.mlsd.side_effect = side_eff_ftp_mlsd

    def mock_retr_binary(command: str, callback: object):
        logging.info(f"Faking the download of {command}")
        try:
            with data_file.open("rb") as data_fh:
                callback(data_fh.read())
        except OSError as err:
            raise ftp_error_reply from err

    mock_ftp.retrbinary.side_effect = mock_retr_binary
    
    with exception:
        connected_ftp = establish_ftp(mock_ftp, ftp_url, accession)
        mock_ftp.connect.assert_called()
        mock_ftp.login.assert_called()
        mock_ftp.cwd.assert_called_with("genomes/all/GCA/017/607/445")
    
    with exception:
        download_files(connected_ftp, accession, data_dir, max_redo)
        mock_ftp.cwd.assert_called()
        mock_ftp.mlsd.assert_called()
        mock_ftp.retrbinary.assert_called()

#################
@pytest.mark.dependency(name="test_get_files_selection")
@pytest.mark.parametrize(
    "download_dir, files_expected, expectation",
    [
        pytest.param(
            False,
            {},
            pytest.raises(FileDownloadError),
            id="Error case, data dir not provided",
        ),
        pytest.param(
            True,
            {
            "report" : "GCA_017607445.1_ASM1760744v1_assembly_report.txt",
            "fasta_dna" : "GCA_017607445.1_ASM1760744v1_genomic.fna.gz",
            "fasta_pep" : "GCA_017607445.1_ASM1760744v1_protein.faa.gz",
            "gff3_raw" : "GCA_017607445.1_ASM1760744v1_genomic.gff.gz",
            "gbff" : "GCA_017607445.1_ASM1760744v1_genomic.gbff.gz",
            },
            does_not_raise(),
            id="Normal case, data dir provided",
        ),
    ],
)
def test_get_files_selection(data_dir: Path, download_dir: bool, files_expected: dict, expectation: ContextManager) -> None:
    """Test the get a subset of downloaded files function `get_files_selection()`
    
    Args:
        download_dir: Path to specific location of downloaded files.
    """

    if download_dir:
        download_dir = data_dir
    else:
        download_dir = Path("")

    with expectation:
        subset_files = get_files_selection(download_dir)
        for file_end_name in subset_files.keys():
            expected_file = files_expected[file_end_name]
            test_data_file_name = subset_files[file_end_name].split("/")[-1]
            assert test_data_file_name == expected_file

##################
# @pytest.mark.dependency(name="test_download_all_files")
@pytest.mark.parametrize(
    "accession, is_dir, files_downloaded, md5_return, exception",
    [
        pytest.param(
            "GCA_017607445.1",
            True,
            {
            "report" : "GCA_017607445.1_ASM1760744v1_assembly_report.txt",
            "fasta_dna" : "GCA_017607445.1_ASM1760744v1_genomic.fna.gz",
            "fasta_pep" : "GCA_017607445.1_ASM1760744v1_protein.faa.gz",
            "gff3_raw" : "GCA_017607445.1_ASM1760744v1_genomic.gff.gz",
            "gbff" : "GCA_017607445.1_ASM1760744v1_genomic.gbff.gz",
            },
            True,
            does_not_raise(),
            id="Case 1: Good accession, dir exists",
        ),
        pytest.param(
            "GCA_017607445.1",
            False,
            {
            "report" : "GCA_017607445.1_ASM1760744v1_assembly_report.txt",
            "fasta_dna" : "GCA_017607445.1_ASM1760744v1_genomic.fna.gz",
            "fasta_pep" : "GCA_017607445.1_ASM1760744v1_protein.faa.gz",
            "gff3_raw" : "GCA_017607445.1_ASM1760744v1_genomic.gff.gz",
            "gbff" : "GCA_017607445.1_ASM1760744v1_genomic.gbff.gz",
            },
            False,
            pytest.raises(FileExistsError),
            id="Case 2: Good accession, Not a dir",
        ),
        pytest.param(
            "GCA_017607445.1",
            True,
            {},
            False,
            pytest.raises(FileDownloadError),
            id="Bad case: Download files list empty",
        ),
    ]
)
@patch("ftplib.FTP")
@patch("os.mkdir")
@patch("ensembl.io.genomio.assembly.download.get_files_selection")
@patch("ensembl.io.genomio.assembly.download.download_files")
@patch("ensembl.io.genomio.assembly.download._download_file")
@patch("ensembl.io.genomio.assembly.download.md5_files")
def test_retrieve_assembly_data(mock_retrieve: Mock, mock_download_singlefile: Mock, 
                                mock_download_files: Mock, mock_file_select: Mock, 
                                mock_os: MagicMock, mock_ftp: Mock, data_dir: Path, 
                                accession: str, is_dir: bool, files_downloaded: dict, 
                                md5_return: bool, exception: ContextManager
                                ) -> None:
    """Test of main assembly downloading function:

    Args:
        accession:
        is_dir:
        files_downloaded:
        expectation:
    """

    if is_dir == False:
        download_dir = Path(data_dir / str("NotADir"))
    else:
        download_dir = data_dir

    def side_eff_conn(url:str):
        if not url:
            raise Exception()

    mock_ftp.connect.side_effect = side_eff_conn
    
    def mock_mkdir(command: str):
        logging.info(f"Faking the creation of directory {command} on path {download_dir}")
        try:
            return download_dir
        except OSError as err:
            raise FileExistsError from err
        
    mock_os.mkdir.side_effect = mock_mkdir
    mock_file_select.return_value = files_downloaded
    mock_retrieve.return_value = md5_return

    with exception:
        retrieve_assembly_data(accession, download_dir)
        print(f"retrieve_assembly_data({accession}, {download_dir})")
        assert mock_os.mkdir.called_with(download_dir)
        assert mock_os.md5_files.called_once_with("md5checksums.txt")
        assert mock_download_files.download_files.called_once()
        assert mock_download_singlefile._download_file.called_once()
        assert mock_file_select.get_files_selection.called_with(files_downloaded)


# Before adding test on file_selection:
# /homes/lcampbell/MyDevelopmentSpace/2024_hackathons/EnsMetazoaHackathon_Jan24/ensembl-genomio/src/python/ensembl/io/genomio/assembly/download.py       172     67    61%   158-187, 221-224, 231, 239, 241-242, 250, 268-273, 289-290, 312-345, 350-359

# After adding test on file_selection:
# /homes/lcampbell/MyDevelopmentSpace/2024_hackathons/EnsMetazoaHackathon_Jan24/ensembl-genomio/src/python/ensembl/io/genomio/assembly/download.py       172     59    66%   158-187, 221-224, 231, 239, 241-242, 250, 312-345, 350-359

# After adding test for could not download file (_download_file)
# /homes/lcampbell/MyDevelopmentSpace/2024_hackathons/EnsMetazoaHackathon_Jan24/ensembl-genomio/src/python/ensembl/io/genomio/assembly/download.py       164     51    69%   158-187, 221-224, 231, 239, 241-242, 312-345

# After adding inital tests on main download_files()
# /homes/lcampbell/MyDevelopmentSpace/2024_hackathons/EnsMetazoaHackathon_Jan24/ensembl-genomio/src/python/ensembl/io/genomio/assembly/download.py       164     26    84%   223, 230, 238, 240-241, 311-344

## Post adding proper test on download_files
# /homes/lcampbell/MyDevelopmentSpace/2024_hackathons/EnsMetazoaHackathon_Jan24/ensembl-genomio/src/python/ensembl/io/genomio/assembly/download.py       164     26    84%   223, 230, 238, 240-241, 311-344

## Latest changes to download_files      
# /homes/lcampbell/MyDevelopmentSpace/2024_hackathons/EnsMetazoaHackathon_Jan24/ensembl-genomio/src/python/ensembl/io/genomio/assembly/download.py       164     40    76%   172-188, 222-225, 232, 240,242-243, 313-346