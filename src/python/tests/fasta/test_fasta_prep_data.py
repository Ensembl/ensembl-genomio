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
"""Unit testing of :mod:`ensembl.io.genomio.fasta.process` module.

"""

import filecmp
from pathlib import Path
from typing import List

import pytest

import ensembl.io.genomio.fasta.process as FastaProcessing


class TestFastaProcess:
    """Tests for the Fasta module, submodule 'process'."""

    tmp_dir: Path
    data_dir: Path

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, tmp_dir: Path, files_dir: Path):
        """Loads necessary fixtures and values as class attributes."""
        type(self).tmp_dir = tmp_dir
        type(self).data_dir = files_dir / "process_fasta"

    @pytest.mark.parametrize(
        "input_fasta, input_gbff, expected_output_fasta, pep_mode",
        [
            (
                "input.protein.fa.gz",
                "input.gbff.gz",
                "output.protein.fa",
                True,
            ),
            (
                "input.protein.fa.gz",
                None,
                "output.protein.fa",
                True,
            ),
            (
                "input.protein.fa.gz",
                None,
                "output.protein.fa",
                False,
            ),
            ("input.fna.gz", None, "output.nuc.fna", False),
            ("input.fna.gz", None, "output.nuc.fna", True),
        ],
    )
    def test_fasta_prep(
        self, tmp_path: Path, input_fasta: str, input_gbff: str, expected_output_fasta: str, pep_mode: bool
    ) -> None:
        """Tests the `process.prep_fasta_data()` function.

        Args:
            tmp_path: Where temporary files will be created.
            input_fasta: Name of the fasta file with example input, in the test folder.
            input_gbff: Name of the input GBFF example input, in the test folder.
            expected_output_fasta: Name of the output fasta file with expected output, in the test folder.
            pep_mode: Boolean flag to set processing in peptide mode.
        """

        fasta_input_path = self.data_dir / input_fasta
        if input_gbff is not None:
            gbff_input_path = self.data_dir / input_gbff
        else:
            gbff_input_path = input_gbff
        fasta_output_path = tmp_path / f"{input_fasta}.test.fasta"

        FastaProcessing.prep_fasta_data(
            fasta_infile=fasta_input_path,
            genbank_infile=gbff_input_path,
            fasta_outfile=fasta_output_path,
            peptide_mode=pep_mode,
        )

        expected_path = self.data_dir / expected_output_fasta

        assert filecmp.cmp(fasta_output_path, expected_path)

    @pytest.mark.parametrize(
        "input_gbff, excluded_seq_region, output, scenario",
        [
            (
                "input.gbff.gz",
                ["LR605957.1"],
                ["VWP78966.1", "VWP78967.1", "VWP78968.1"],
                "Normal",
            ),
            (  
                "input.mod.gbff.gz", 
                ["LR605957.1"], 
                ["VWP78966.1", "VWP78967.1", "VWP78968.1"],
                "Modified"
            ),
        ],
    )
    def test_exclude_seq_regions(
        self,
        input_gbff: str,
        excluded_seq_region: List[str],
        output: List[str],
        scenario: str,
    ) -> None:
        """Tests the `process.get_peptides_to_exclude()` function.

        Args:
            input_gbff: Name of the input GBFF example input, in the test folder.
            excluded_seq_region: Set of sequence regions to be excluded.
            output: Set of proteins expected to be excluded given excluded_seq_region
            scenario: Normal or modified function test.
        """

        if input_gbff is not None:
            gbff_input_path = self.data_dir / input_gbff
        else:
            gbff_input_path = input_gbff

        seq_region_to_exclude = set(excluded_seq_region)

        # Testing 'Normal' operation i.e. intact and correct GenBank gbff
        if scenario == "Normal":
            excluded_proteins = FastaProcessing.get_peptides_to_exclude(
                gbff_input_path, seq_region_to_exclude
            )
            assert set(excluded_proteins) == set(output)  # Expecting 3 peptides on seq_region (LR605957.1)
        # Improper operation with corrupt GenBank gbff
        else:
            pytest.raises(
                FastaProcessing.FastaParserError,
                FastaProcessing.get_peptides_to_exclude,
                gbff_input_path,
                seq_region_to_exclude,
            )
