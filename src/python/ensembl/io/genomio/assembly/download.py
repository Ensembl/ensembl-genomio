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
"""Download an assembly data files from INSDC or RefSeq."""

__all__ = [
    "FileDownloadError",
    "FTPConnectionError",
    "UnsupportedFormatError",
    "md5_files",
    "get_checksums",
    "download_files",
    "get_files_selection",
    "get_root_name",
    "retrieve_assembly_data",
]

from ftplib import FTP
import hashlib
from importlib import reload
import logging
from os import PathLike
from pathlib import Path
import re
import time
from typing import Dict

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

_FILE_ENDS = {
    "assembly_report.txt": "report",
    "genomic.fna.gz": "fasta_dna",
    "protein.faa.gz": "fasta_pep",
    "genomic.gff.gz": "gff3_raw",
    "genomic.gbff.gz": "gbff",
}


class FileDownloadError(Exception):
    """When a file download fails or there is a problem with that file."""

class FTPConnectionError(Exception):
    """Error while initialising an FTP connection."""

class UnsupportedFormatError(Exception):
    """When a string does not have the expected format."""

class FTPConnection:
    """Custom Class for establishing an FTP connection."""

    def establish_ftp (ftp_conn: FTP, ftp_url: str, sub_dir_path: Path) -> FTP:
        """Function to establish an FTP connection (FTP Class).

            Args:
                ftp_conn: FTP Class object
                ftp_url: specific FTP url in connection request
                sub_dir_path: Path of sub directory housing required data for download.

            Returns:
                An open connection via FTP() class object
        """

        try:
            ftp_conn.connect(ftp_url)
            ftp_conn.login()
            ftp_conn.cwd(str(sub_dir_path))
            return ftp_conn
        except:
            raise FTPConnectionError(f"Could not create FTP connection on {ftp_url} remote path {sub_dir_path}")


def md5_files(dl_dir: Path, md5: Path) -> bool:
    """
    Check all files checksums with the sums listed in a checksum file, if available.
    Return False if there is no checksum file, or a file is missing, or has a wrong checksum.
    """
    # Get or set md5 file to user or default setting
    if md5 is None:
        md5_file = "md5checksums.txt"
    else:
        md5_file = md5
        
    # Get checksums and compare
    md5_path = dl_dir / md5_file
    sums = get_checksums(md5_path)
    if not sums:
        return False
    logging.info(f" File sums from {md5_path}: {len(sums)}")
    for dl_file, checksum in sums.items():
        for end in _FILE_ENDS:
            if dl_file.endswith(end) and not dl_file.endswith(f"_from_{end}"):
                file_path = dl_dir / dl_file
                if not file_path.is_file():
                    logging.warning(f" No file {file_path} found")
                    return False
                # Check the file checksum
                with file_path.open(mode="rb") as f:
                    content = f.read()
                    file_sum = hashlib.md5(content).hexdigest()
                if file_sum != checksum:
                    logging.warning(f" File {file_path} checksum doesn't match")
                    return False
                logging.info(f" File checksum ok {file_path}")
    logging.info(" All checksums OK")
    return True


def get_checksums(checksum_path: Path) -> Dict[str, str]:
    """
    Get a dict of checksums from a file, with file names as keys and sums as values
    """
    sums: Dict[str, str] = {}
    if not checksum_path.is_file():
        return sums
    with checksum_path.open(mode="r") as fh:
        for line in fh:
            checksum, file_path = line.strip().split("  ")
            file_path = file_path[2:]
            if not file_path.find("/") >= 0:
                sums[file_path] = checksum
    return sums

def download_files(accession: str, dl_dir: Path, max_redo: int) -> None:
    """
    Given an INSDC accession, download all available files from the ftp to the download dir
    """
    match = re.match(r"(GC[AF])_([0-9]{3})([0-9]{3})([0-9]{3})\.?([0-9]+)", accession)
    if not match:
        raise UnsupportedFormatError(f"Could not recognize GCA accession format: {accession}")
    gca = match.group(1)
    part1 = match.group(2)
    part2 = match.group(3)
    part3 = match.group(4)

    # Establish connection to remote FTP server
    ftp_url = "ftp.ncbi.nlm.nih.gov"
    sub_dir = Path("genomes", "all", gca, part1, part2, part3)
    ftp_instance = FTP()
    ftp_connection = FTPConnection.establish_ftp(ftp_instance, ftp_url, sub_dir)

    # # Get the list of assemblies for this accession
    for ftp_dir, _ in ftp_connection.mlsd():
        if re.match(accession, ftp_dir):
            ftp_connection.cwd(ftp_dir)
            # First, get the md5sum file
            md5_file = "md5checksums.txt"
            md5_path = dl_dir / md5_file
            with md5_path.open("wb") as fp:
                ftp_connection.retrbinary(f"RETR {md5_file}", fp.write)
            md5_sums = get_checksums(md5_path)

            # Get all the files
            for ftp_file, _ in ftp_connection.mlsd():
                for end in _FILE_ENDS:
                    if ftp_file.endswith(end) and not ftp_file.endswith(f"_from_{end}"):
                        _download_file(ftp_connection, ftp_file, md5_sums, dl_dir, max_redo)


def _download_file(
    ftp_conn: FTP, ftp_file: str, md5_sums: Dict[str, str], dl_dir: Path, max_redo: int = 0
) -> None:
    """TODO"""
    has_md5 = True
    expected_sum = ""
    if not ftp_file in md5_sums:
        logging.warning(f" File not in the md5 checksums: {ftp_file}")
        has_md5 = False
    else:
        expected_sum = md5_sums[ftp_file]
    local_path = Path(dl_dir, ftp_file)

    # File exists? Check md5sum before anything else
    if local_path.is_file():
        if has_md5:
            with local_path.open(mode="rb") as fp:
                content = fp.read()
                file_sum = hashlib.md5(content).hexdigest()
                if file_sum == expected_sum:
                    logging.info(f" File {local_path} is already downloaded properly")
                    return
        else:
            logging.info(f" Can't check file (no md5sum), using it as is: {local_path}")
    file_sum = ""
    redo = 0

    while (file_sum != expected_sum) and (redo <= max_redo):
        redo += 1
        if redo > 1:
            time.sleep(3)

        # Download the file
        logging.info(f" Downloading file {ftp_file}, try {redo}...")
        try:
            with local_path.open(mode="wb") as fp:
                ftp_conn.retrbinary(f"RETR {ftp_file}", fp.write)
        except EOFError:
            continue
        if not has_md5:
            file_sum = ""
            continue
        # Compute checksum
        with local_path.open(mode="rb") as fp:
            content = fp.read()
            file_sum = hashlib.md5(content).hexdigest()
    if expected_sum == file_sum:
        logging.info(f" Downloaded file properly to {local_path}")
    else:
        raise FileDownloadError(f"Could not download file {ftp_file} after {redo} tries")


def get_files_selection(dl_dir: Path) -> Dict[str, str]:
    """
    Among all the files downloaded, only keep a subset for which we use a controlled name.
    Return a dict[name] = file_path
    The file_path is relative to the download dir
    Current names are defined in _FILE_ENDS
    """
    files = {}
    root_name = get_root_name(dl_dir)
    if root_name == "":
        raise FileDownloadError(f"Could not determine the files root name in {dl_dir}")
    for dl_file in dl_dir.iterdir():
        for end, name in _FILE_ENDS.items():
            file_with_end = dl_file.name.endswith(end) and not dl_file.name.endswith(f"_from_{end}")
            if (root_name and dl_file.name == root_name + end) or file_with_end:
                files[name] = str(dl_file)
    return files


def get_root_name(dl_dir: Path) -> str:
    """Get root name for assembly files, using the report file as base"""
    root_name = ""
    for dl_file in dl_dir.iterdir():
        matches = re.search("^(.+_)assembly_report.txt", dl_file.name)
        if matches:
            root_name = matches.group(1)
            break
    return root_name


def retrieve_assembly_data(
    accession: str,
    download_dir: PathLike,
    max_increment: int = 0,
    max_redo: int = 3,
) -> None:
    """TODO

    Args:
        accession: Genome Assembly accession
        download_dir: Path to directory used to store retrieved
        max_increment: If you want to allow assembly versions
        max_redo: Set max number of times to retry downloading a file

    """
    download_dir = Path(download_dir)

    # Configure logging
    log_file = f"{accession}_download.log"
    reload(logging)
    logging.basicConfig(
        filename=log_file, format="%(levelname)s:%(message)s", filemode="w", level=logging.DEBUG
    )

    # Set and create dedicated dir for download
    if not download_dir.is_dir():
        download_dir.mkdir(parents=True)

    # Download if files don't exist or fail checksum
    if not md5_files(download_dir, None):
        logging.info(" Download the files")

        for increment in range(0, max_increment + 1):
            if increment > 0:
                logging.info(f" Increment accession version once from {accession}")
                version = int(accession[-1])
                version += 1
                accession = accession[:-1] + str(version)
                download_dir.mkdir(parents=True, exist_ok=True)
            download_files(accession, download_dir, max_redo)

        if not md5_files(download_dir, None):
            raise FileDownloadError("Failed md5sum of downloaded files")

    # Select specific files and give them a name
    files = get_files_selection(download_dir)

    if len(files) == 0:
        raise FileDownloadError("No file downloaded")


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description="Download an assembly data files from INSDC or RefSeq.")
    parser.add_argument("--accession", required=True, help="Genome assembly accession")
    parser.add_argument_dst_path(
        "--download_dir", default=Path.cwd(), help="Folder where the data will be downloaded"
    )
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    retrieve_assembly_data(args.accession, args.download_dir)
