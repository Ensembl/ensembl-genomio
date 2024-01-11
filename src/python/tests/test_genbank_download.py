"""Unit testing of 'ensembl.io.genomio.genbank.download' module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_genbank_download.py

"""
import filecmp
from pathlib import Path
import pytest
import requests

from ensembl.io.genomio.genbank.download import download_genbank

class TestGenbank_download:
    """Tests for the 'download_genbank' class"""   
        
    data_dir: Path
    tmp_dir: Path

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, tmp_dir: Path, files_dir: Path):
        """Loads necessary fixtures and values as class attributes."""
        type(self).tmp_dir = tmp_dir
        type(self).data_dir = files_dir /"genbank_download"

    def test_successful_download(self):
        """Tests the 'download_genbank()' method.
    
        Args:
        accession: Genbank accession to be downloaded
        output: path to the expected output file 
        """

        accession = "CM023531.1"
        output_file = self.tmp_dir/f"{accession}.gb"
        download_genbank(accession, output_file)
        expected_output = self.data_dir/f"{accession}.gb"
        assert filecmp.cmp(output_file, expected_output)
