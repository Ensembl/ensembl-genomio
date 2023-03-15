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
import os
import requests


class DownloadError(Exception):
    """In case a download failed."""

    def __init__(self, msg):
        self.msg = msg


class download_genbank(eHive.BaseRunnable):
    def param_defaults(self):
        return {}

    def run(self):
        accession = self.param_required("gb_accession")
        main_download_dir = self.param("download_dir")
        download_dir = main_download_dir + "/" + accession

        # Set and create dedicated dir for download
        if not os.path.isdir(download_dir):
            os.makedirs(download_dir)

        # Download the file
        gb_path = self.download_genbank(accession, download_dir)

        output = {"gb_file": gb_path, "gb_accession": accession}
        self.dataflow(output, 2)

    @staticmethod
    def download_genbank(accession: str, dl_dir: str) -> str:
        """
        Given a GenBank accession, download the corresponding file in GenBank format
        """
        dl_path = os.path.join(dl_dir, accession + ".gb")

        # Don't redownload the file
        if os.path.exists(dl_path):
            return dl_path

        # Get the list of assemblies for this accession
        e_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        e_params = {
            "db": "nuccore",
            "rettype": "gb",
            "retmode": "text",
        }
        e_params["id"] = accession

        result = requests.get(e_url, params=e_params)

        if result and result.status_code == 200:
            with open(dl_path, "wb") as gff:
                gff.write(result.content)
            print(f"GFF file write to {dl_path}")
            return dl_path
        else:
            raise DownloadError(f"Could not download the genbank file: {result}")
