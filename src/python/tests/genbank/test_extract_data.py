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
import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation,SimpleLocation
from unittest.mock import MagicMock, patch
from Bio import SeqIO
from BCBio import GFF

from ensembl.io.genomio.genbank.extract_data import (
    FormattedFilesGenerator, GBParseError
)

class TestWriteFormattedFiles:

    @pytest.fixture
    def formatted_files_generator(self,data_dir: Path):
        """Call the function `FormattedFilesGenerator` with set parameters"""
        prod_name= "TEST_prod"
        gb_file = "input_file.gb"
        gb_file_path = data_dir/gb_file
        prefix = "TEST"
        return FormattedFilesGenerator(prod_name, gb_file_path, prefix)
    
    def test_parse_genbank(self, data_dir: Path, formatted_files_generator) -> None:
        """Test that parse_genbank method correctly parses genbank files."""
        gb_file_path = data_dir / formatted_files_generator.gb_file
        formatted_files_generator.parse_genbank(gb_file_path)
        assert len(formatted_files_generator.seq_records) >= 1

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

    def test_format_genome_json_with_production_name(self, formatted_files_generator, tmp_path: Path):
        """Test write_genome_json generates correct json output."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(FeatureLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]})        
        record1.features.append(gene_feature1)

        record2 = SeqRecord(Seq("ATGC"), id="record2")
        gene_feature1 = SeqFeature(FeatureLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrB"]})        
        record1.features.append(gene_feature1)
    
        formatted_files_generator.seq_records = [record1, record2]
        formatted_files_generator.files["genome"] = tmp_path / "genome.json"

        formatted_files_generator._format_genome_data()
        assert (tmp_path / "genome.json").exists()
        with open(tmp_path / "genome.json", "r") as f:
            genome_data = json.load(f)
            assert genome_data["species"]["production_name"] == "TEST_prod"
            assert genome_data["assembly"]["accession"] == "GCA_000000000"
            assert genome_data["assembly"]["version"] == 1
            assert genome_data["added_seq"]["region_name"] == ["record1", "record2"]
        
    def test_write_fasta_dna(self, tmp_path: Path, formatted_files_generator) -> None:
        """Check fasta DNA sequences are written as expected."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(FeatureLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]})        
        record1.features.append(gene_feature1)
        formatted_files_generator.seq_records = [record1]

        formatted_files_generator.files["fasta_dna"] = tmp_path / "test.fasta"

        formatted_files_generator._write_fasta_dna()

        assert (tmp_path / "test.fasta").exists()
        fasta_pep= SeqIO.read((tmp_path / "test.fasta"), "fasta")
        assert fasta_pep.id == record1.id
        assert fasta_pep.seq == record1.seq

    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._get_codon_table")
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._prepare_location")
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._write_seq_region_json")
    def test_format_seq_region_json(self, mock_codon_table, mock_org_location, mock_write_seq_json, formatted_files_generator) -> None:
        """Check _seq_region_json is written as expected."""
        record1 = SeqRecord(Seq("ATGC"), id="record1", annotations={"topology":"circular"})
        CDS_feature1 = SeqFeature(FeatureLocation(10,20), type="CDS", qualifiers= {"gene": ["GlyrA"], "transl_table": "2"}) 
        record1.features.append(CDS_feature1)
        record1.organelle = "mitochondrion"
        formatted_files_generator.seq_records = [record1]
        
        formatted_files_generator._format_write_seq_json()
        
        # Call the method to be tested
        mock_codon_table.assert_called_once()
        mock_org_location.assert_called_once()
        mock_write_seq_json.assert_called_once()

    def test_format_write_genes_gff(self, formatted_files_generator 
    ):
        """Check gene features in GFF3 format are generated as expected."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(FeatureLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]})        
        record1.features.append(gene_feature1)
        
        formatted_files_generator.seq_records = [record1]
        mock_all_ids = ["ID1", "ID2", "ID3"]
        mock_new_record = record1
        mock_peptides = ["pep1", "pep2"]
        
        # Set the mock data to be returned by _parse_record and mock all the methods
        formatted_files_generator._parse_record = MagicMock(return_value=(
            mock_new_record,  # Mock the new record
            mock_all_ids,     # Mock the list of all IDs
            mock_peptides     # Mock the list of peptides
        ))    
        formatted_files_generator._write_genes_gff = MagicMock()
        formatted_files_generator._write_pep_fasta = MagicMock()

        # Call the method to be tested
        formatted_files_generator._format_write_genes_gff()
        
        # Assert that _parse_record was called
        formatted_files_generator._parse_record.assert_called_once()
        formatted_files_generator._write_genes_gff.assert_called_once()
        formatted_files_generator._write_pep_fasta.assert_called_once()

        mock_all_ids_with_duplicates = ["ID1", "ID2", "ID1"]
        formatted_files_generator._parse_record = MagicMock(return_value=(
            mock_new_record,  # Mock the new record
            mock_all_ids_with_duplicates,     # Mock the list of all IDs
            mock_peptides     # Mock the list of peptides
        ))  
        with pytest.raises(GBParseError):
            formatted_files_generator._format_write_genes_gff()

    def test_write_genes_gff(self, formatted_files_generator, tmp_path: Path):
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(FeatureLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]}) 
        CDS_feature1 = SeqFeature(FeatureLocation(10,15), type="CDS", qualifiers= {"gene": ["GlyrA"], "transl_table": "2"}) 
        record1.features=[gene_feature1,CDS_feature1]
        formatted_files_generator.seq_records = [record1]

        formatted_files_generator.files["gene_models"] = tmp_path / "genes.gff"
        
        formatted_files_generator._write_genes_gff(formatted_files_generator.seq_records)
        assert (formatted_files_generator.files["gene_models"]).exists()
        for rec in GFF.parse(formatted_files_generator.files["gene_models"]):
            assert rec.id == record1.id
            assert len(rec.features)== len(record1.features)

    def test_write_pep_fasta(self, formatted_files_generator, tmp_path: Path):
        record1 = SeqRecord(Seq("MFLRTQARFFHATTKKM"), id="cds-record1")
        CDS_feature1 = SeqFeature(FeatureLocation(10,20), type="CDS", qualifiers= {"gene": ["GlyrA"], "transl_table": "2"}) 
        record1.features.append(CDS_feature1)
        formatted_files_generator.files["fasta_pep"] = tmp_path / "pep.fasta"
        
        formatted_files_generator._write_pep_fasta(record1)
        assert (tmp_path / "pep.fasta").exists()

        fasta_pep= SeqIO.read((tmp_path / "pep.fasta"), "fasta")
        assert fasta_pep.id == record1.id
        assert fasta_pep.seq == record1.seq