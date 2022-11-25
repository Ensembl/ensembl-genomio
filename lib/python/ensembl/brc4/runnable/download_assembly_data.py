#!env python3
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

from pathlib import Path
from typing import Any, Dict
import eHive
import os
import re
import ftplib
import hashlib
from pathlib import Path

FILE_ENDS = {
    "assembly_report.txt": "report",
    "genomic.fna.gz": "fasta_dna",
    "protein.faa.gz": "fasta_pep",
    "genomic.gff.gz": "gff3_raw",
    "genomic.gbff.gz": "gbff",
}
    

class download_assembly_data(eHive.BaseRunnable):

    @staticmethod
    def param_defaults() -> Dict[str, Any]:
        return {
            # Set this manually to a higher value if you want to allow assembly versions
            # higher than the one provided (it will fail if the version given is not the latest)
            "max_increment": 0,
            # Set max number of times to retry downloading a file
            "max_redo": 3,
        }

    def run(self):
        accession = self.param_required('accession')
        main_download_dir = Path(self.param_required('download_dir'))
        download_dir = main_download_dir / accession

        # Set and create dedicated dir for download
        if not download_dir.is_dir():
            os.makedirs(download_dir)

        # Download if files don't exist or fail checksum
        if not self.md5_files(download_dir):
            print("Download the files")
            
            max_increment = self.param('max_increment')

            for increment in range(0, max_increment + 1):
                if increment > 0:
                    print(f"Increment accession version once from {accession}")
                    version = int(accession[-1])
                    version += 1
                    accession = accession[:-1] + str(version)
                    download_dir = main_download_dir / accession
                    if not download_dir.is_dir():
                        download_dir.mkdir(parents=True)
                self.download_files(accession, download_dir)

            if not self.md5_files(download_dir):
                raise Exception("Failed md5sum of downloaded files")

        # Select specific files and give them a name
        files = self.get_files_selection(download_dir)

        if len(files) == 0:
            raise Exception("No file downloaded")

        # Output all those named files + dir
        self.dataflow(files, 2)

    def md5_files(self, dl_dir: Path) -> bool:
        """
        Check all files checksums with the sums listed in a checksum file, if available.
        Return False if there is no checksum file, or a file is missing, or has a wrong checksum.
        """
        md5_file = "md5checksums.txt"

        # Get checksums and compare
        md5_path = dl_dir / md5_file
        sums = self.get_checksums(md5_path)
        
        if not sums:
            return

        print(f"File sums from {md5_path}: {len(sums)}")
        
        for dl_file, checksum in sums.items():
            for end in FILE_ENDS:
                if dl_file.endswith(end) and not dl_file.endswith(f"_from_{end}"):
                    file_path = dl_dir / dl_file
                    
                    if not file_path.is_file():
                        print(f"No file {file_path} found")
                        return False
                    
                    # Check the file checksum
                    with file_path.open(mode='rb') as f:
                        content = f.read()
                        file_sum = hashlib.md5(content).hexdigest()
                    if file_sum != checksum:
                        print(f"File {file_path} checksum doesn't match")
                        return False
                    else:
                        print(f"File checksum ok {file_path}")

        print("All checksums OK")
        return True
    
    @staticmethod
    def get_checksums(checksum_path: Path) -> Dict[str, str]:
        """
        Get a dict of checksums from a file, with file names as keys and sums as values
        """

        if not checksum_path.is_file():
            return
        
        sums = {}
        with checksum_path.open(mode='r') as fh:
            for line in fh:
                checksum, file_path = line.strip().split("  ")
                file_path = file_path[2:]
                if not file_path.find("/") >= 0:
                    sums[file_path] = checksum

        return sums

    def download_files(self, accession: str, dl_dir: Path) -> None:
        """
        Given an INSDC accession, download all available files from the ftp to the download dir
        """
        match = re.match(r'(GC[AF])_([0-9]{3})([0-9]{3})([0-9]{3})\.?([0-9]+)', accession)
        gca = match.group(1)
        part1 = match.group(2)
        part2 = match.group(3)
        part3 = match.group(4)
        parts = (gca, part1, part2, part3)

        # Get the list of assemblies for this accession
        ftp_url = "ftp.ncbi.nlm.nih.gov"
        sub_dir = Path("genomes", "all", gca, part1, part2, part3)
        f = ftplib.FTP()
        f.connect(ftp_url)
        f.login()
        f.cwd(str(sub_dir))

        max_redo = self.param("max_redo")

        for (ftp_dir, entry) in f.mlsd():
            if re.match(accession, ftp_dir):
                f.cwd(ftp_dir)
            
                # First, get the md5sum file
                md5_file = 'md5checksums.txt'
                md5_path = dl_dir / md5_file
                with md5_path.open("wb") as fp:
                    f.retrbinary(f"RETR {md5_file}", fp.write)
                md5_sums = self.get_checksums(md5_path)
                
                # Get all the files
                for (ftp_file, file_entry) in f.mlsd():
                    has_md5 = True
                    expected_sum = ''
                    for end in FILE_ENDS:
                        if ftp_file.endswith(end) and not ftp_file.endswith(f"_from_{end}"):
                            if not ftp_file in md5_sums:
                                print(f"File not in {md5_file}: {ftp_file}")
                                has_md5 = False
                            else:
                                expected_sum = md5_sums[ftp_file]
                            local_path = Path(dl_dir, ftp_file)

                            # File exists? Check md5sum before anything else
                            if local_path.is_file():
                                if has_md5:
                                    with local_path.open(mode='rb') as fp:
                                        content = fp.read()
                                        file_sum = hashlib.md5(content).hexdigest()
                                        if file_sum == expected_sum:
                                            print(f"File {local_path} is already downloaded properly")
                                            continue
                                else:
                                    print(f"Can't check file (no md5sum), using it as is: {local_path}")

                            file_sum = ''
                            redo = 0
                            while (file_sum != expected_sum) and (redo <= max_redo):
                                redo += 1
                                print(f"Downloading file {ftp_file}, try {redo}...")

                                # Download the file
                                try:
                                    with local_path.open(mode='wb') as fp:
                                        f.retrbinary(f"RETR {ftp_file}", fp.write)
                                except EOFError:
                                    continue
                                
                                if not has_md5:
                                    file_sum = ''
                                    continue

                                # Compute checksum
                                with local_path.open(mode='rb') as fp:
                                    content = fp.read()
                                    file_sum = hashlib.md5(content).hexdigest()
                            if expected_sum == file_sum:
                                print(f"Downloaded file properly to {local_path}")
                            else:
                                raise Exception(f"Could not download file {ftp_file} after {redo} tries")


    def get_files_selection(self, dl_dir: Path) -> Dict[str, str]:
        """
        Among all the files downloaded, only keep a subset for which we use a controlled name.
        Return a dict[name] = file_path
        The file_path is relative to the download dir
        
        Current names are defined in FILE_ENDS
        """
        files = {}

        root_name = self.get_root_name(dl_dir)
        if root_name == '':
            raise Exception(f"Could not determine the files root name in {dl_dir}")

        for dl_file in dl_dir.iterdir():
            for end, name in FILE_ENDS.items():
                file_with_end = dl_file.name.endswith(end) and not dl_file.name.endswith(f"_from_{end}")
                if (root_name and dl_file.name == root_name + end) or file_with_end:
                    files[name] = str(dl_file)
        return files

    @staticmethod
    def get_root_name(dl_dir: Path) -> str:
        """Get root name for assembly files, using the report file as base"""

        root_name = ''
        for dl_file in dl_dir.iterdir():
            matches = re.search("^(.+_)assembly_report.txt", dl_file.name)
            if matches:
                root_name = matches.group(1)
                break
        
        return root_name
