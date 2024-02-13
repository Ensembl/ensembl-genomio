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
"""Unit testing of `ensembl.io.genomio.genbank.extract_data` module.

Typical usage example::
    $ pytest test_extract_data_files.py

"""
from pathlib import Path
import filecmp
import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation

from ensembl.io.genomio.genbank.extract_data import (
    FormattedFilesGenerator, GBParseError
)

class TestWriteFormattedFiles:

    @pytest.fixture
    def formatted_files_generator(self,data_dir):
        """Call the function `FormattedFilesGenerator` with set parameters"""
        prod_name= "TEST_prod"
        gb_file = "input_file.gb"
        gb_file_path = data_dir/gb_file
        prefix = "TEST"
        return FormattedFilesGenerator(prod_name, gb_file_path, prefix)
    
    @pytest.mark.parametrize(
        "expected",
        [
            ({"S43128.1": "mitochondrion"}),
        ],
    )
    def test_get_organella(self, data_dir: Path, expected: str, formatted_files_generator) -> None:
        """Test that organellas are correctly identified."""
        gb_file_path = data_dir / formatted_files_generator.gb_file
        organella = formatted_files_generator._get_organella(gb_file_path)
        assert organella == expected

    def test_parse_genbank(self, data_dir: Path, formatted_files_generator) -> None:
        """Test that parse_genbank method correctly parses genbank files."""
        gb_file_path = data_dir / formatted_files_generator.gb_file
        formatted_files_generator.parse_genbank(gb_file_path)
        assert len(formatted_files_generator.seq_records) >= 1
     
    def test_write_genome_json(self, data_dir: Path, tmp_path: Path, formatted_files_generator) -> None:
        """Test write_genome_json generates correct json output."""
        out_dir = tmp_path
        formatted_files_generator.extract_gb(out_dir)
        output= tmp_path/ "genome.json"
        expected_path = data_dir/ f"output_genome.json"
        assert filecmp.cmp(output, expected_path)

    def test_write_fasta_dna(self, data_dir: Path, tmp_path: Path, formatted_files_generator) -> None:
        """Check fasta DNA sequences are written as expected."""
        out_dir = tmp_path
        formatted_files_generator.extract_gb(out_dir)
        output= tmp_path/ "dna.fasta"
        expected_path = data_dir/ f"output_dna.fasta"
        assert filecmp.cmp(output, expected_path)

    def test_write_seq_region_json(self, data_dir: Path, tmp_path: Path, formatted_files_generator) -> None:
        """Check _seq_region_json is written as expected."""
        out_dir = tmp_path
        formatted_files_generator.extract_gb(out_dir)
        output= tmp_path/ "seq_region.json"
        expected_path = data_dir/ f"output_seq_region.json"
        assert filecmp.cmp(output, expected_path)
    
    def test_write_genes_gff(self, assert_files, data_dir: Path, tmp_path: Path,  formatted_files_generator, ):
        """Check gene features in GFF3 format are generated as expected."""
        peptides=[]
        formatted_files_generator.files["gene_models"] = tmp_path / "genes.gff"
        output= formatted_files_generator.files["gene_models"]
        expected_path = data_dir/ f"output_genes.gff"
        if peptides:
            formatted_files_generator.files["fasta_pep"] = tmp_path / "pep.fasta"
            output_pep =  formatted_files_generator.files["fasta_pep"] 
            expected_pep = data_dir/ f"output_pep.fasta"
            assert_files(output_pep, expected_pep)

        assert_files(output, expected_path)

    def test_duplicate_features(self, formatted_files_generator, tmp_path):
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene")
        record1.features.append(gene_feature1)

        record2 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature2 = SeqFeature(SimpleLocation(10,30), type="gene")
        record2.features.append(gene_feature2)

        formatted_files_generator.seq_records = [record1, record2]
        formatted_files_generator.files["gene_models"] = tmp_path / "genes.gff"
        
        with pytest.raises(GBParseError):
            formatted_files_generator._write_genes_gff()
    



