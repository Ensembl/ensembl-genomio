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
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
import filecmp
import pytest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


from ensembl.io.genomio.genbank.extract_data import (
    FormattedFilesGenerator, UnsupportedData
)

class TestFormattedFilesGenerator:

    @pytest.fixture
    def formatted_files_generator(data_dir):
        prod_name= "TEST_prod"
        gb_file = "input_file.gb"
        gb_file_path = data_dir/gb_file
        prefix = "TEST"
        return FormattedFilesGenerator(prod_name, gb_file_path, prefix)

    def test_parse_genbank(data_dir, formatted_files_generator):
        gb_file_path = data_dir / f"input_file.gb"
        formatted_files_generator.parse_genbank(gb_file_path)
        assert len(formatted_files_generator.seq_records) >= 1
     
    def test_write_genome_json(data_dir, tmp_path, formatted_files_generator):
        """Test that organellas are correctly identified."""
        out_dir = tmp_path
        formatted_files_generator.extract_gb(out_dir)
        output= tmp_path/ "genome.json"
        expected_path = data_dir/ f"output_genome.json"
        assert filecmp.cmp(output, expected_path)

    def test_write_fasta_dna(data_dir, tmp_path, formatted_files_generator):
        """Test that organellas are correctly identified."""
        out_dir = tmp_path
        formatted_files_generator.extract_gb(out_dir)
        output= tmp_path/ "dna.fasta"
        expected_path = data_dir/ f"output_dna.fasta"
        assert filecmp.cmp(output, expected_path)

    def test_write_seq_region_json(tmp_path, data_dir, formatted_files_generator):
        out_dir = tmp_path
        formatted_files_generator.extract_gb(out_dir)
        output= tmp_path/ "seq_region.json"
        expected_path = data_dir/ f"output_seq_region.json"
        assert filecmp.cmp(output, expected_path)

