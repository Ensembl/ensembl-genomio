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
import json
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

    def test_write_genome_json_with_production_name(self, formatted_files_generator, tmp_path, caplog):
        """Test write_genome_json generates correct json output."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]})        
        record1.features.append(gene_feature1)

        record2 = SeqRecord(Seq("ATGC"), id="record2")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrB"]})        
        record1.features.append(gene_feature1)
    
        formatted_files_generator.seq_records = [record1, record2]
        formatted_files_generator.files["genome"] = tmp_path / "genome.json"

        formatted_files_generator._write_genome_json()

        assert (tmp_path / "genome.json").exists()
        with open(tmp_path / "genome.json", "r") as f:
            genome_data = json.load(f)
            assert genome_data["species"]["production_name"] == "TEST_prod"
            assert genome_data["assembly"]["accession"] == "GCA_000000000"
            assert genome_data["assembly"]["version"] == 1
            assert genome_data["added_seq"]["region_name"] == ["record1", "record2"]
        
    def test_write_fasta_dna(self, data_dir: Path, tmp_path: Path, formatted_files_generator) -> None:
        """Check fasta DNA sequences are written as expected."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]})        
        record1.features.append(gene_feature1)
        formatted_files_generator.seq_records = [record1]

        formatted_files_generator.files["fasta_dna"] = tmp_path / "test.fasta"

        formatted_files_generator._write_fasta_dna()

        assert (tmp_path / "test.fasta").exists()

    def test_write_seq_region_json(self, data_dir: Path, tmp_path: Path, formatted_files_generator) -> None:
        """Check _seq_region_json is written as expected."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"], "transl_table": "2"})        
        record1.features.append(gene_feature1)
        formatted_files_generator.seq_records = [record1]
        formatted_files_generator.files["seq_region"] = tmp_path / "seq_region.json"
        formatted_files_generator._write_seq_region_json()
        assert (tmp_path / "seq_region.json").exists()
        with open(tmp_path / "seq_region.json", "r") as f:
            seq_region = json.load(f)
            assert len(seq_region) == 1
            assert seq_region["codon_table"] == "2"
    
    
    def test_write_genes_gff(self, assert_files, data_dir: Path, tmp_path: Path,  formatted_files_generator, ):
        """Check gene features in GFF3 format are generated as expected."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]})        
        record1.features.append(gene_feature1)
    
        formatted_files_generator.seq_records = [record1]

        formatted_files_generator.files["gene_models"] = tmp_path / "genes.gff"
        formatted_files_generator._write_genes_gff()
        
        assert (tmp_path / "genes.gff").exists()

    def test_duplicate_features(self, formatted_files_generator, tmp_path):
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene")
        record1.features.append(gene_feature1)

        record2 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature2 = SeqFeature(SimpleLocation(10,30), type="rna")
        record2.features.append(gene_feature2)

        formatted_files_generator.seq_records = [record1, record2]
        formatted_files_generator.files["gene_models"] = tmp_path / "genes.gff"
        
        with pytest.raises(GBParseError):
            formatted_files_generator._write_genes_gff()
    



