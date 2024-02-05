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
"""Unit testing of `ensembl.io.genomio.genome_stats.compare` module.

Typical usage example::
    $ pytest test_compare.py

"""

import filecmp
import pytest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


from ensembl.io.genomio.genbank.extract_data import (
    FormattedFilesGenerator,
)

class TestFormattedFilesGenerator:

    prod_name= "TEST_prod"
    gb_file = "input_file.gb"

    @pytest.mark.parametrize("expected",
                        [
                            ({"NC_000884.1": "mitochondrion"}),
                        ],
                        )
    def test_get_organella(self,data_dir, expected):
        """Test that organellas are correctly identified."""

        gb_file_path = data_dir / self.gb_file
        gb = FormattedFilesGenerator(self.prod_name,self.gb_file)
        organella = gb._get_organella(gb_file_path)
        assert organella == expected

    def test_parse_genbank(self,data_dir):
        gb_file_path = data_dir / self.gb_file
        gb = FormattedFilesGenerator(self.prod_name,self.gb_file)
        gb.parse_genbank(gb_file_path)
        assert len(gb.seq_records) >= 1

    def test_write_genome_json(self, data_dir, tmp_path):
        """Test that organellas are correctly identified."""

        gb_file_path = data_dir / self.gb_file
        out_dir = tmp_path
        gb = FormattedFilesGenerator(self.prod_name,gb_file_path)
        gb.extract_gb(out_dir)
        output= tmp_path/ "genome.json"
        expected_path = data_dir/ f"output_genome.json"
        assert filecmp.cmp(output, expected_path)

    def test_write_fasta_dna(self, data_dir, tmp_path):
        """Test that organellas are correctly identified."""

        gb_file_path = data_dir / self.gb_file
        out_dir = tmp_path
        gb = FormattedFilesGenerator(self.prod_name,gb_file_path)
        gb.extract_gb(out_dir)
        output= tmp_path/ "dna.fasta"
        expected_path = data_dir/ f"output_dna.fasta"
        assert filecmp.cmp(output, expected_path)

    def test_get_codon_table(self, data_dir):
        gb_file_path = data_dir / self.gb_file
        gb = FormattedFilesGenerator(self.prod_name,self.gb_file)
        gb.parse_genbank(gb_file_path)
        for seq in gb.seq_records:
            codon_table = gb._get_codon_table(seq)
        assert codon_table == '2'

    def test_write_seq_region_json(self, tmp_path, data_dir):
        gb_file_path = data_dir / self.gb_file
        out_dir = tmp_path
        gb = FormattedFilesGenerator(self.prod_name,gb_file_path)
        gb.extract_gb(out_dir)
        output= tmp_path/ "seq_region.json"
        expected_path = data_dir/ f"output_seq_region.json"
        assert filecmp.cmp(output, expected_path)

