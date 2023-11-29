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
"""Unit testing of :mod:`ensembl.io.genomio.gff3.process_gff3` module.

"""

from typing import Dict, Union
import pytest
from unittest.mock import MagicMock

from ensembl.io.genomio.database import CoreServer


TEST_LIST = [
    "speciesA_genus_core_60_110_1",
    "speciesB_genus_sp_core_60_110_1",
    "prefix_speciesC_genus_core_60_110_1",
    "prefix_speciesD_genus_sp_core_60_110_1",
    "speciesE_genus_core_61_110_1",
    "speciesF_genus_core_61_111_1",
]


class TestCoreServer:
    """Tests for the database.core_server module."""

    @pytest.mark.parametrize(
        "parameters, output_size",
        [
            ({}, len(TEST_LIST)),
            ({"build": 60}, 4),
            ({"build": 61}, 2),
            ({"prefix": "foobar"}, 0),
            ({"prefix": "prefix_"}, 2),
            ({"version": 110}, 5),
            ({"version": 111}, 1),
            ({"version": 999}, 0),
            ({"dbname_re": r""}, len(TEST_LIST)),
            ({"dbname_re": r"foobar"}, 0),
            ({"dbname_re": r"^species._.*"}, 4),
            ({"dbname_re": r"species._.*_sp_"}, 2),
            ({"dbname_re": r"^prefix_"}, 2),
            ({"prefix": "prefix_", "build": 60}, 2),
            ({"build": 61, "version": 110}, 1),
        ],
    )
    def test_get_cores(self, parameters: Dict[str, Union[str, int]], output_size: int):
        """Test the CoreServer.get_cores() method."""

        # Fake server with mock get_all_core_names()
        server_url = "sqlite:///:memory:"
        server = CoreServer(server_url)
        server.get_all_core_names = MagicMock(return_value=TEST_LIST)

        # NB: get_cores() gets its data from the mocked get_all_core_names()
        all_cores = server.get_cores(**parameters)
        assert len(all_cores) == output_size
