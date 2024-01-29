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

from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from contextlib import nullcontext as does_not_raise
from typing import ContextManager

import pytest
from ftplib import FTP


from ensembl.io.genomio.assembly.download import FTPConnection, FTPConnectionError


class TestDownloadAssembly:
    """Tests class for `download_assembly` module."""

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
        self,
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
