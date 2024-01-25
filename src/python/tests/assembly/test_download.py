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
    $ pytest test_download.py

"""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from pytest import raises

from ensembl.io.genomio.assembly.download import download_files, md5_files

@pytest.mark.parametrize(
    "accession",
    [
        ("GCA_000002765.3"),
    ],
)
@patch("ensembl.io.genomio.genbank.download.requests.get")