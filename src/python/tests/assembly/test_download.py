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
"""Unit testing of `ensembl.io.genomio.assembly.download` module."""
# pylint: disable=too-many-positional-arguments

from contextlib import nullcontext as does_not_raise
import filecmp
import logging
from pathlib import Path
from typing import Callable, ContextManager, Optional
from unittest.mock import Mock, patch, MagicMock

from ftplib import error_reply as ftp_error_reply
import pytest

from ensembl.io.genomio.assembly.download import UnsupportedFormatError
from ensembl.io.genomio.assembly.download import FileDownloadError
from ensembl.io.genomio.assembly.download import FTPConnectionError
from ensembl.io.genomio.assembly.download import establish_ftp
from ensembl.io.genomio.assembly.download import get_checksums, md5_files, _download_file
from ensembl.io.genomio.assembly.download import download_files, get_files_selection, retrieve_assembly_data


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
            "ftp.ncbi.nlm.nih.gov",
            "gcx_017607445.1",
            pytest.raises(UnsupportedFormatError),
            id="Failed connection improper INSDC accession",
        ),
    ],
)
@patch("ensembl.io.genomio.assembly.download.FTP", autospec=True)
def test_ftp_connection(
    mock_ftp: Mock,
    ftp_url: str,
    accession: str,
    expectation: ContextManager,
) -> None:
    """Tests the FTPConnection method `establish_ftp()`.

    Args:
        mock_ftp: Mock FTP object.
        ftp_url: FTP URL.
        sub_dir: Subdirectory path.
        expectation: Context manager expected raise exception.
    """

    def side_eff_conn(url: str) -> None:
        if not url:
            raise FTPConnectionError()

    mock_ftp.connect.side_effect = side_eff_conn

    with expectation:
        connected_ftp = establish_ftp(mock_ftp, ftp_url, accession)
        connected_ftp.connect.assert_called_once_with(ftp_url)  # type: ignore[attr-defined]
        connected_ftp.login.assert_called_once()  # type: ignore[attr-defined]
        connected_ftp.cwd.assert_called_once_with("genomes/all/GCA/017/607/445")  # type: ignore[attr-defined]


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
def test_checksums(
    data_dir: Path, checksum_file: Path, checksum: Optional[str], expectation: ContextManager
) -> None:
    """Tests the `download.get_checksums()` function.

    Args:
        data_dir: Path to test data root dir
        checksum_file: File name containing checksums
        checksum: Test MD5 checksum
        expectation: Context manager expected raise exception
    """
    with expectation:
        md5_input_path = data_dir / checksum_file
        obtained_checksums = get_checksums(md5_input_path)
        assert obtained_checksums.get("annotation_hashes.txt", None) == checksum


#################
@pytest.mark.dependency(name="test_md5_files", depends=["test_checksums"])
@pytest.mark.parametrize(
    "md5_file, md5_path, checksum_bool",
    [
        pytest.param("md5checksums.txt", None, True, id="Normal case"),
        pytest.param("wrong_md5_checksums.txt", None, False, id="Incorrect md5 checksum"),
        pytest.param(None, None, True, id="No md5file specified, resort to default"),
        pytest.param(None, Path("*"), False, id="Incompatible os path '*'"),
        pytest.param("missing_file_md5.txt", None, False, id="md5 checksum with ref of missing file"),
    ],
)
def test_md5_files(data_dir: Path, md5_file: str, md5_path: Optional[Path], checksum_bool: bool) -> None:
    """Tests the md5_files() function
    Args:
        data_dir: Path to test data root dir
        md5_file: MD5 file used for test
        md5_path: Path location to MD5 file
        checksum_bool: Test comparison checksum value
    """
    if md5_file is None:
        return_bool_on_md5files = md5_files(data_dir, md5_path=md5_path)
    else:
        return_bool_on_md5files = md5_files(data_dir, md5_path=md5_path, md5_filename=md5_file)
    assert return_bool_on_md5files == checksum_bool


#################
@pytest.mark.dependency(name="test_download_single_file")
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
        pytest.param(
            "test_ftp_file.txt",
            {"test_ftp_file.txt": "e98b980b442fdb2a21877dcc55e11848"},
            does_not_raise(),
            id="Max retry attempts non-zero",
        ),
        pytest.param(
            "test_ftp_file.txt",
            {"test_ftp_file.txt": "e98b980b442fdb2a21877XxxXxxXxxXx"},
            pytest.raises(FileDownloadError),
            id="Incorrect md5sum",
        ),
    ],
)
@patch("ensembl.io.genomio.assembly.download.FTP", autospec=True)
def test_download_single_file(
    mock_ftp: Mock,
    data_dir: Path,
    tmp_path: Path,
    ftp_file: str,
    md5_sums: dict,
    expectation: ContextManager,
) -> None:
    """Tests the private function _download_file.

    Args:
        mock_ftp: Mock FTP object
        data_dir: Path to test data root dir
        tmp_path: Temp dir created for test
        ftp_file: FTP file which to mock download
        md5_sums: FTP file and md5_sum value pair
        expectation: Context manager expected raise exception
    """

    data_file = data_dir / ftp_file
    retr_file = tmp_path / ftp_file

    def mock_retr_binary(command: str, callback: Callable) -> None:
        logging.info(f"Faking the download of {command}")
        try:
            with data_file.open("rb") as data_fh:
                callback(data_fh.read())
        except OSError as err:
            raise ftp_error_reply from err

    mock_ftp.retrbinary.side_effect = mock_retr_binary

    with expectation:
        _download_file(mock_ftp, ftp_file, md5_sums, tmp_path)
        assert filecmp.cmp(data_file, retr_file)


#################
@pytest.mark.dependency(name="test_download_all_files", depends=["test_download_single_file"])
@pytest.mark.parametrize(
    "ftp_url, ftp_accession, compare_accession, md5, exception, max_redo",
    [
        pytest.param(
            "ftp.ncbi.nlm.nih.gov",
            "GCA_017607445.1",
            "GCA_017607445.1",
            "md5checksums_ftp_source.txt",
            does_not_raise(),
            2,
            id="Normal case: Properly formatted GC[A/F] accession",
        ),
        pytest.param(
            "ftp.ncbi.nlm.nih.gov",
            "GCA_017607445.1",
            "GCA_544706710.1",
            "md5checksums_ftp_source.txt",
            does_not_raise(),
            2,
            id="Error case: GCA doesn't match ftp_dir path",
        ),
    ],
)
@patch("ensembl.io.genomio.assembly.download.FTP", autospec=True)
def test_download_all_files(
    mock_ftp: MagicMock,
    data_dir: Path,
    ftp_url: str,
    ftp_accession: str,
    compare_accession: str,
    md5: str,
    exception: ContextManager,
    max_redo: int,
) -> None:
    """Tests the download.download_files() function

    Args:
        mock_ftp: Mock of `ensembl.io.genomio.assembly.download.FTP` object
        data_dir: Path to local test data folder
        ftp_url: Test param to specify ftp connection url
        ftp_accession: Defines expected accession
        compare_accession: Defines test of expected accession
        md5: Source file for md5 checksums to inspect
        expectation: Context manager expected raise exception
    """

    data_file = data_dir / md5

    def side_eff_ftp_mlsd() -> list[tuple[str, list[str]]]:
        mlsd_ret = []
        files = data_dir.glob("*GCA_017607445.1*.gz")
        for file_path in files:
            ftp_file = Path(file_path).name
            mlsd_ret.append((str(ftp_file), ["fact_1", "fact_2"]))

        return mlsd_ret

    mock_ftp.mlsd.side_effect = side_eff_ftp_mlsd

    def mock_retr_binary(command: str, callback: Callable) -> None:
        logging.info(f"Faking the download of {command}")
        try:
            with data_file.open("rb") as data_fh:
                callback(data_fh.read())
        except OSError as err:
            raise ftp_error_reply from err

    mock_ftp.retrbinary.side_effect = mock_retr_binary

    with exception:
        connected_ftp = establish_ftp(mock_ftp, ftp_url, ftp_accession)
        mock_ftp.connect.assert_called()
        mock_ftp.login.assert_called()
        mock_ftp.cwd.assert_called_with("genomes/all/GCA/017/607/445")

        download_files(connected_ftp, compare_accession, data_dir, max_redo)
        mock_ftp.cwd.assert_called()
        mock_ftp.mlsd.assert_called()

        if ftp_accession == compare_accession:
            mock_ftp.retrbinary.assert_called()
        else:
            mock_ftp.retrbinary.assert_not_called()


#################
@pytest.mark.dependency(name="test_get_files_selection")
@pytest.mark.parametrize(
    "has_download_dir, files_expected, expectation",
    [
        pytest.param(
            True,
            {
                "report": "GCA_017607445.1_ASM1760744v1_assembly_report.txt",
                "fasta_dna": "GCA_017607445.1_ASM1760744v1_genomic.fna.gz",
                "fasta_pep": "GCA_017607445.1_ASM1760744v1_protein.faa.gz",
                "gff3_raw": "GCA_017607445.1_ASM1760744v1_genomic.gff.gz",
                "gbff": "GCA_017607445.1_ASM1760744v1_genomic.gbff.gz",
            },
            does_not_raise(),
            id="Normal case, data dir provided",
        ),
        pytest.param(
            False,
            {},
            pytest.raises(FileDownloadError),
            id="Error case, data dir not provided",
        ),
    ],
)
def test_get_files_selection(
    data_dir: Path, has_download_dir: bool, files_expected: dict, expectation: ContextManager
) -> None:
    """Tests the `download.get_files_selection()` function.

    Args:
        download_dir: Path to specific location of downloaded files.
        files_expected: Defines contents of test files downloaded
        expectation: Context manager expected raise exception
    """

    if has_download_dir:
        download_dir = data_dir
    else:
        download_dir = Path()

    with expectation:
        subset_files = get_files_selection(download_dir)
        for file_end_name, file_path in subset_files.items():
            expected_file = files_expected[file_end_name]
            test_data_file_name = Path(file_path).name
            assert test_data_file_name == expected_file


##################
@pytest.mark.dependency(name="test_retrieve_assembly_data")
@pytest.mark.parametrize(
    "accession, is_dir, files_downloaded, md5_return, exception",
    [
        pytest.param(
            "GCA_017607445.1",
            True,
            {
                "report": "GCA_017607445.1_ASM1760744v1_assembly_report.txt",
                "fasta_dna": "GCA_017607445.1_ASM1760744v1_genomic.fna.gz",
                "fasta_pep": "GCA_017607445.1_ASM1760744v1_protein.faa.gz",
                "gff3_raw": "GCA_017607445.1_ASM1760744v1_genomic.gff.gz",
                "gbff": "GCA_017607445.1_ASM1760744v1_genomic.gbff.gz",
            },
            True,
            does_not_raise(),
            id="Case 1: Good accession, dir exists",
        ),
        pytest.param(
            "GCA_017607445.1",
            False,
            0,
            False,
            pytest.raises(FileExistsError),
            id="Case 2: Good accession, Not a dir",
        ),
        pytest.param(
            "GCA_017607445.1",
            True,
            0,
            False,
            pytest.raises(FileDownloadError),
            id="Case 3: Download files list empty",
        ),
    ],
)
@patch("ensembl.io.genomio.assembly.download.FTP", autospec=True)
@patch("ensembl.io.genomio.assembly.download.get_files_selection")
@patch("ensembl.io.genomio.assembly.download.download_files")
@patch("ensembl.io.genomio.assembly.download._download_file")
@patch("ensembl.io.genomio.assembly.download.md5_files")
def test_retrieve_assembly_data(
    mock_retrieve: Mock,
    mock_download_single_file: Mock,
    mock_download_files: Mock,
    mock_file_select: Mock,
    mock_ftp: Mock,
    data_dir: Path,
    accession: str,
    is_dir: bool,
    files_downloaded: dict,
    md5_return: bool,
    exception: ContextManager,
) -> None:
    """Test of master function download.retrieve_assembly_data() which calls sub
    functions for downloading assembly data files.

    Args:
        accession: Accession of desired genome assembly
        is_dir: Param to define state of result output dir
        files_downloaded: Defines contents of test files marked as downloaded
        expectation: Context manager expected raise exception
    """

    if is_dir:
        download_dir = data_dir
    else:
        download_dir = data_dir / "test_ftp_file.txt"

    def side_eff_conn(url: str) -> None:
        if not url:
            raise FileDownloadError()

    mock_ftp.connect.side_effect = side_eff_conn
    mock_file_select.return_value = files_downloaded
    mock_retrieve.return_value = md5_return

    with exception:
        retrieve_assembly_data(accession, download_dir, 2)
        mock_file_select.assert_called_with(download_dir)
        # Download is never called with the current examples
        mock_download_files.assert_not_called()
        mock_download_single_file.assert_not_called()
        # mock_download_files.assert_called_once()
        # mock_download_single_file.assert_called_once()
