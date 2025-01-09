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
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse
from io import StringIO
import pytest
from unittest.mock import MagicMock, patch
from pytest import TempPathFactory

from ensembl.utils.archive import open_gz_file

from ensembl.io.genomio.fasta.compare import CompareFasta
from ensembl.io.genomio.fasta.compare import SeqGroup

class TestWriteFormattedFiles:
    """Test if all the expected output files are generated and formatted correctly"""

    fasta1_path = "fasta1"
    fasta2_path = "fasta2"
    output_path = "output"

    @pytest.fixture(scope="class", autouse=True)
    def compare_fasta_instance(self) -> CompareFasta:
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
    
    @patch("ensembl.io.genomio.fasta.compare.SeqIO.parse")
    @patch("ensembl.io.genomio.fasta.compare.open_gz_file")
    def test_read_fasta(self, mock_open_gz_file, mock_seqio_parse, compare_fasta_instance):
        """
        Tests the `read_fasta` method of CompareFasta.
    
        This test ensures that `read_fasta` correctly parses a FASTA file,
        replaces non-CGTA characters with 'N', and returns the expected dictionary.
        """
    
        # Define the mock FASTA file path
        fasta_path = Path("dummy1.fasta")
    
        # Create mock SeqRecord objects
        seq_record1 = SeqRecord(Seq("ATCGNNNN"), id="seq1")
        seq_record2 = SeqRecord(Seq("NANTGCGT"), id="seq2")
    
        # Expected output after replacing non-CGTA characters with 'N'
        expected_output = {
            "seq1": "ATCGNNNN",
            "seq2": "NANTGCGT",
        }
    
        # Configure the mock for open_gz_file to return a mock file handle
        mock_fasta_fh = StringIO(">seq1\nATCGNNNN\n>seq2\nNANTGCGT\n")
        mock_open_gz_file.return_value = mock_fasta_fh
    
        # Configure the mock for SeqIO.parse to return the mock SeqRecord objects
        mock_seqio_parse.return_value = [seq_record1, seq_record2]
    
        # Call the method under test
        result = compare_fasta_instance.read_fasta(fasta_path)
    
        # Assertions to ensure methods were called correctly
        mock_open_gz_file.assert_called_once_with(fasta_path)
        mock_seqio_parse.assert_called_once_with(mock_fasta_fh, "fasta")
    
        # Assert that the result matches the expected output
        assert result == expected_output, "The read_fasta method did not return the expected dictionary."

    @pytest.mark.parametrize(
    "seq_dict1, seq_dict2, expected_common, expected_comp",
    [
        pytest.param(
            # Input dictionaries with SeqGroup objects
            {"AATCG": SeqGroup("id1")},
            {"AATCG": SeqGroup("id2")},
            # Expected common output
            {"id1": "id2"},
            # Expected comparison logs
            [],
            id="Matching single-ID groups"
        ),
        pytest.param(
            # Multi-ID groups with identical counts
            {"ATTCG": SeqGroup("id1"), "ATCG": SeqGroup("id4"), "ATCG": SeqGroup("id5")},
            {"ATCG": SeqGroup("id2"), "ATTCG": SeqGroup("id3")},
            {'id1': 'id3', 'id5': 'id2 OR id3', 'id4': 'id2 OR id3'},
            ["Matched 2 identical groups of sequences: id5, id4 and id2, id3"],
            id="Matching multi-ID groups with same count"
        ),
    ]
    )
    def test_find_common_groups(self,seq_dict1, seq_dict2, expected_common, expected_comp, compare_fasta_instance) -> None:
        """Tests the `find_common_groups` function."""

        # Add extra IDs for multi-ID cases to mimic realistic scenarios
        if "ATCG" in seq_dict1 :
            seq_dict1["ATCG"].add_id("id4")
        if "ATCG" in seq_dict2 :
            seq_dict2["ATCG"].add_id("id3")

        common, comp = compare_fasta_instance.find_common_groups(seq_dict1,seq_dict2)
        assert common == expected_common
        assert comp == expected_comp

    @pytest.mark.parametrize(
    "only1, only2, expected_comp",
    [
        pytest.param(
            # Input dictionaries with SeqGroup objects
            {"AATCGN": "id1"},
            {"AATCGNN": "id2"},
            # Expected comparison logs
            ["Please check extra Ns added in your fasta2 in id1 and id2"],
            id="Matching single-ID groups with extra Ns"
        ),
        pytest.param(
            # Input dictionaries with SeqGroup objects
            {"AATCGN": "id1"},
            {"AATCGG": "id2"},
            # Expected comparison logs
            ["ALERT INSERTIONS at the end or diff assembly level id1 and id2",
            "sequences have the same length, check id1 and id2"],
            id="Matching single-ID groups with different sequences"
        ),
    ]
    )
    def test_compare_seq_for_Ns(self, only1, only2, expected_comp, compare_fasta_instance):
        """Tests the `compare_seq_for_Ns` function."""

         # Clear the comp list to isolate this test
        compare_fasta_instance.comp = []
        compare_fasta_instance.compare_seq_for_Ns(only1, only2)
        assert compare_fasta_instance.comp == expected_comp

    @pytest.mark.parametrize(
        "mock_seq1, mock_seq2, mock_seq1_dict, mock_seq2_dict",
        [
        pytest.param(
            ["ACTG", "CGTA"],
            ["ACTG", "CGTANN"],
            {"ACTG": "id1", "CGTA": "id2"},
            {"ACTG": "id3", "CGTANN": "id4"},
            id="Case with differing sequence lengths"
        ),
    ],
    )
    @patch("ensembl.io.genomio.fasta.compare.CompareFasta.write_results")
    @patch("ensembl.io.genomio.fasta.compare.CompareFasta.build_seq_dict")
    @patch("ensembl.io.genomio.fasta.compare.CompareFasta.read_fasta")
    def test_compare_seqs(self,
        mock_read_fasta,
        mock_build_seq_dict,
        mock_write_results,
        compare_fasta_instance,
        mock_seq1,
        mock_seq2,
        mock_seq1_dict,
        mock_seq2_dict,
    ):
        """Tests the `compare_seqs` function."""

        # Configure the mocks
        mock_read_fasta.side_effect = [mock_seq1, mock_seq2]
        mock_build_seq_dict.side_effect = [mock_seq1_dict, mock_seq2_dict]

        # Call the method
        compare_fasta_instance.compare_seqs()

        # Assertions
        mock_read_fasta.assert_called()
        mock_build_seq_dict.assert_called()
        mock_write_results.assert_called()

        # Check the comparison logs
        assert len(compare_fasta_instance.comp) > 0

    def test_write_results(self, compare_fasta_instance, tmp_path_factory: TempPathFactory):
        """Tests the `write_results` method."""
        
        temp = tmp_path_factory.mktemp("temp")
        compare_fasta_instance.output_dir = temp
        compare_fasta_instance.write_results()
        expected_file = Path(temp) / "compare.log"  # Update filename as per your implementation
        assert expected_file.exists(), f"Expected file {expected_file} does not exist."
    
        