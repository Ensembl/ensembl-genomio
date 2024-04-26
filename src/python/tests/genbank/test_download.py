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
"""Unit testing of `ensembl.io.genomio.genbank.download` module.

The unit testing is divided into one test class per submodule/class found in this module, and one test method
per public function/class method.

Typical usage example::
    $ pytest test_download.py

"""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from pytest import raises

from ensembl.io.genomio.genbank.download import download_genbank, DownloadError


@pytest.mark.parametrize(
    "accession",
    [
        ("CM023531.1"),
    ],
)
@patch("ensembl.io.genomio.genbank.download.requests.get")
class TestDownloadGenbank:
    """Tests for the `download_genbank` class"""

    def test_successful_download(self, mock_requests_get: Mock, tmp_path: Path, accession: str) -> None:
        """Tests the successful download of `download_genbank()` method.

        Args:
            mock_requests_get: A mock of `request.get()` method.
            tmp_path: Function-scoped temporary directory fixture.
            accession: Genbank accession to be downloaded.
        """

        # Set success_code and content as an attribute to the mock object
        mock_requests_get.return_value.status_code = 200
        mock_content = b"The genbank download for the following accession"
        mock_requests_get.return_value.content = mock_content

        # Temporary location where we want to store the mock output file
        output_file = tmp_path / f"{accession}.gb"
        download_genbank(accession, output_file)

        # Checking if the url has been called
        mock_requests_get.assert_called_once_with(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={"db": "nuccore", "rettype": "gbwithparts", "retmode": "text", "id": accession},
            timeout=60,
        )

        # Assert that the content was written to the temporary file
        with open(output_file, "rb") as f:
            file_content = f.read()
        assert file_content == mock_content

    def test_failed_download(self, mock_requests_failed: Mock, tmp_path: Path, accession: str) -> None:
        """Tests the downloading failure of `download_genbank()` method.

        Args:
            mock_requests_failed: A mock of `request.get()` method.
            tmp_path: Function-scoped temporary directory fixture.
            accession: Genbank accession to be downloaded.
        """

        output_file = tmp_path / f"{accession}.gb"
        # Set the mock status code to 404 for request not found
        mock_requests_failed.return_value.status_code = 404
        # Raise an error
        with raises(DownloadError):
            download_genbank(accession, output_file)
