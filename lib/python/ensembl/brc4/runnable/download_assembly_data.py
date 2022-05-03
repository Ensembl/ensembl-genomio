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



import eHive
import json
import os, re
import ftplib
import hashlib

class download_assembly_data(eHive.BaseRunnable):
    
    @staticmethod
    def param_defaults():
        return {
            # Set this manually to a higher value if you want to allow assembly versions
            # higher than the one provided (it will fail if the version given is not the latest)
            "max_increment" : 0,
        }

    def run(self):
        accession = self.param_required('accession')
        main_download_dir = self.param('download_dir')
        download_dir = main_download_dir + "/" + accession

        # Set and create dedicated dir for download
        if not os.path.isdir(download_dir):
            os.makedirs(download_dir)

        # Download if files don't exist or fail checksum
        if not self.md5_files(download_dir):
            print("Download the files")
            
            max_increment = self.param('max_increment')

            for increment in range(0, max_increment + 1):
                if increment > 0:
                    print("Increment accession version once from %s" % accession)
                    version = int(accession[-1])
                    version += 1
                    accession = accession[:-1] + str(version)
                    download_dir = main_download_dir + "/" + accession
                    if not os.path.isdir(download_dir):
                        os.makedirs(download_dir)
                try:
                    self.download_files(accession, download_dir)
                    break
                except Exception:
                    print("Can't download files for %s" % accession)

            if not self.md5_files(download_dir):
                raise Exception("Failed md5sum of downloaded files")

        # Select specific files and give them a name
        files = self.get_files_selection(download_dir)

        if len(files) == 0:
            raise Exception("No file downloaded")

        # Output all those named files + dir
        self.dataflow(files, 2)

    @staticmethod
    def get_json(json_path):
        """
        Retrieve the json data from a json file
        """
        with open(json_path) as json_file:
            json_data = json.load(json_file)
        return json_data

    def md5_files(self, dl_dir):
        """
        Check all files checksums with the sums listed in a checksum file, if available.
        Return False if there is no checksum file, or a file is missing, or has a wrong checksum.
        """
        md5_file = "md5checksums.txt"

        # Get checksums and compare
        md5_path = os.path.join(dl_dir, md5_file)
        sums = self.get_checksums(md5_path)
        
        if not sums:
            return

        print("File sums from %s: %d" % (md5_path, len(sums)))
        
        for dl_file, checksum in sums.items():
            file_path = os.path.join(dl_dir, dl_file)
            
            if not os.path.isfile(file_path):
                print("No file %s found" % file_path)
                return False
            
            # Check the file checksum
            with open(file_path, mode='rb') as f:
                content = f.read()
                file_sum = hashlib.md5(content).hexdigest()
            if file_sum != checksum:
                print("File %s checksum doesn't match" % file_path)
                return False
            else:
                print("File checksum ok %s" % file_path)

        print("All checksums OK")
        return True
    
    @staticmethod
    def get_checksums(checksum_path):
        """
        Get a dict of checksums from a file, with file names as keys and sums as values
        """

        if not os.path.isfile(checksum_path):
            return
        
        sums = {}
        with open(checksum_path, mode='r') as fh:
            for line in fh:
                checksum, file_path = line.strip().split("  ")
                file_path = file_path[2:]
                if not file_path.find("/") >= 0:
                    sums[file_path] = checksum

        return sums

    @staticmethod
    def download_files(accession, dl_dir):
        """
        Given an INSDC accession, download all available files from the ftp to the download dir
        """
        match = re.match("(GC[AF])_([0-9]{3})([0-9]{3})([0-9]{3})\.?([0-9]+)", accession)
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
                    if re.match("^\.", ftp_file):
                        continue
                    if file_entry["type"] == "dir":
                        continue

                    # Copy the file locally
                    local_path = os.path.join(dl_dir, ftp_file)
                    with open(local_path, 'wb') as fp:
                        f.retrbinary("RETR " + ftp_file, fp.write)

    def get_files_selection(self, dl_dir):
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
        if root_name is None:
            raise Exception("Could not determine the files root name in %s" % dl_dir)

        for dl_file in os.listdir(dl_dir):
            for end, name in file_ends.items():
                file_with_end = dl_file.endswith(end) and not dl_file.endswith("_from_" + end)
                if root_name and dl_file == root_name + end or file_with_end:
                    files[name] = os.path.join(dl_dir, dl_file)
        return files

    @staticmethod
    def get_root_name(dl_dir):
        """Get root name for assembly files, using the report file as base"""

        for dl_file in os.listdir(dl_dir):
            matches = re.search("^(.+_)assembly_report.txt", dl_file)
            if matches:
                return matches.group(1)

