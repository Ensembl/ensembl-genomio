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

from unittest.mock import MagicMock

from ensembl.io.genomio.database import CoreServer


TEST_LIST = [
    "speciesA_genus_core_60_110_1",
    "speciesB_genus_sp_core_60_110_1",
    "prefix_speciesC_genus_core_60_110_1",
    "prefix_speciesD_genus_sp_core_60_110_1",
    "speciesE_genus_core_61_110_1",
    "speciesF_genus_sp_core_61_110_1",
]


class TestCoreServer:
    """Tests for the database.core_server module."""

    def test_get_cores(self):
        """Test the CoreServer.get_cores() method."""
        server_url = "sqlite:///:memory:"
        server = CoreServer(server_url)
        server.get_all_core_names = MagicMock(return_value=TEST_LIST)

        all_cores = server.get_cores()
        assert len(all_cores) == len(TEST_LIST)

        prefixed = server.get_cores(prefix="prefix_")
        assert len(prefixed) == 2
