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


class download_assembly_data(eHive.BaseRunnable):
    
    @staticmethod
    def param_defaults() -> Dict[str, Any]:
        return {
            # Set this manually to a higher value if you want to allow assembly versions
            # higher than the one provided (it will fail if the version given is not the latest)
            "max_increment": 0,
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
                    if not os.path.isdir(download_dir):
                        os.makedirs(download_dir)
                try:
                    self.download_files(accession, download_dir)
                    break
                except Exception:
                    print(f"Can't download files for {accession}")

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

        print(f"File sums from {md5_path}:, {len(sums)}")
        
        for dl_file, checksum in sums.items():
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

    @staticmethod
    def download_files(accession: str, dl_dir: Path) -> None:
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
        f = ftplib.FTP()
        f.connect("ftp.ncbi.nlm.nih.gov")
        f.login()
        f.cwd("genomes/all/%s/%s/%s/%s" % parts)

        for (ftp_dir, entry) in f.mlsd():
            if re.match(accession, ftp_dir):
                f.cwd(ftp_dir)
                
                # Get all the files
                for (ftp_file, file_entry) in f.mlsd():
                    if re.match(r'^\.', ftp_file):
                        continue
                    if file_entry["type"] == "dir":
                        continue

                    # Copy the file locally
                    local_path = dl_dir / ftp_file
                    with local_path.open('wb') as fp:
                        f.retrbinary("RETR " + ftp_file, fp.write)

    def get_files_selection(self, dl_dir: Path) -> Dict[str, str]:
        """
        Among all the files downloaded, only keep a subset for which we use a controlled name.
        Return a dict[name] = file_path
        The file_path is relative to the download dir
        
        Current names:
            report
            fasta_dna
            fasta_pep
            gff3
            gbff
        """
        files = {}
        file_ends = {
                "assembly_report.txt" : "report",
                "genomic.fna.gz" : "fasta_dna",
                "protein.faa.gz" : "fasta_pep",
                "genomic.gff.gz" : "gff3_raw",
                "genomic.gbff.gz" : "gbff",
        }

        root_name = self.get_root_name(dl_dir)
        if root_name == '':
            raise Exception(f"Could not determine the files root name in {dl_dir}")

        for dl_file in os.listdir(dl_dir):
            for end, name in file_ends.items():
                file_with_end = dl_file.endswith(end) and not dl_file.endswith("_from_" + end)
                if root_name and dl_file == root_name + end or file_with_end:
                    files[name] = os.path.join(dl_dir, dl_file)
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
