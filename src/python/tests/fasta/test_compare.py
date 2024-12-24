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
"""Unit testing of `ensembl.io.genomio.fasta.compare` module."""
# pylint: disable=too-many-positional-arguments

import pytest

from ensembl.io.genomio.fasta.compare import CompareFasta

class TestWriteFormattedFiles:
    """Test if all the expected output files are generated and formatted correctly"""

    fasta1_path = "fasta1"
    fasta2_path = "fasta2"
    output_path = "output"

    @pytest.fixture(scope="class", autouse=True)
    def compare_fasta_instance(self):
        """Fixture to provide a CompareFasta instance for testing."""

        return CompareFasta(self.fasta1_path, self.fasta2_path, self.output_path)

    @pytest.mark.parametrize(
        "seq_dict, expected",
        [
            pytest.param(
                {"seq1": "ATCG", "seq2": "ATCG"},
                {"ATCG": ["seq1", "seq2"]},
                id="id: seq -> ids added as value when seq match"
            ),
            pytest.param(
                {"seq1": "ATCG", None: "ATATCG"},
                {"ATCG": ["seq1"], "ATATCG": ["None"]},
                id="id: seq -> ids added as value including None"
            ),
            pytest.param(
                {"seq1": "ATTCG", "seq2": "AATCG"},
                {"ATTCG": ["seq1"], "AATCG": ["seq2"]},
                id="id: seq -> ids added as value when no sequences match"
            ),
        ],
    )
    def test_build_seq_dict(self, compare_fasta_instance, seq_dict: dict, expected: dict) -> None:
        """Tests the `build_seq_dict` function.

        Args:
            seq_dict: Input dictionary of sequence IDs to sequences.
            expected: Expected dictionary of sequences mapped to `SeqGroup` objects.
        """
        # Instantiate the CompareFasta object (use dummy paths as required)

        result = compare_fasta_instance.build_seq_dict(seq_dict)

        # Convert `SeqGroup` objects to simple lists for comparison
        result_simple = {seq : seq_ids.ids for seq, seq_ids in result.items()}

        assert result_simple == expected

    




    