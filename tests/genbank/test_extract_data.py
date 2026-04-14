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
"""Unit testing of `ensembl.io.genomio.genbank.extract_data` module."""
# pylint: disable=too-many-positional-arguments

from pathlib import Path
from typing import List
from unittest.mock import Mock, patch

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pytest
from pytest import TempPathFactory

from ensembl.io.genomio.genbank.extract_data import FormattedFilesGenerator, GBParseError
from ensembl.io.genomio.utils import get_json


class TestWriteFormattedFiles:
    """Test if all the expected output files are generated and formatted correctly"""

    prod_name = "TEST_prod"
    gb_file = "input_file.gb"
    prefix = "TEST"

    @pytest.fixture(scope="class", autouse=True)
    def formatted_files_generator(
        self, data_dir: Path, tmp_path_factory: TempPathFactory
    ) -> FormattedFilesGenerator:
        """Call the function `FormattedFilesGenerator` with set parameters.
        Fixture that returns the class of the module that we are testing
        """
        gb_file = self.gb_file
        gb_file_path = data_dir / gb_file
        temp = tmp_path_factory.mktemp("temp")
        return FormattedFilesGenerator(self.prod_name, gb_file_path, self.prefix, out_dir=temp)

    @pytest.mark.dependency(name="parse_genbank")
    def test_parse_genbank(
        self,
        data_dir: Path,
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Test that `parse_genbank()` method to parses the SeqRecords correctly."""
        gb_file_path = data_dir / self.gb_file
        formatted_files_generator.parse_genbank(gb_file_path)
        assert len(formatted_files_generator.seq_records) >= 1

    @pytest.mark.dependency(depends=["parse_genbank"])
    @pytest.mark.parametrize(
        "expected",
        [
            ({"S43128.1": "mitochondrion"}),
        ],
    )
    def test_get_organella(
        self,
        data_dir: Path,
        expected: str,
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Test that organellas are correctly identified in the genbank file"""
        gb_file_path = data_dir / self.gb_file
        # pylint: disable=protected-access
        organella = formatted_files_generator._get_organella(gb_file_path)
        assert organella == expected

    @pytest.mark.dependency(depends=["parse_genbank"])
    def test_write_fasta_dna(
        self,
        tmp_path: Path,
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Check FASTA DNA sequences are written as expected."""
        record = SeqRecord(Seq("ATGC"), id="record")
        gene_feature = SeqFeature(FeatureLocation(10, 20), type="gene", qualifiers={"gene": ["GlyrA"]})
        record.features.append(gene_feature)
        formatted_files_generator.seq_records = [record]

        formatted_files_generator.files["fasta_dna"] = tmp_path / "test.fasta"
        # pylint: disable=protected-access
        formatted_files_generator._write_fasta_dna()

        assert (tmp_path / "test.fasta").exists()
        fasta_pep = SeqIO.read((tmp_path / "test.fasta"), "fasta")
        assert fasta_pep.id == record.id
        assert fasta_pep.seq == record.seq

    @pytest.mark.dependency(depends=["parse_genbank"])
    def test_format_genome_json(
        self,
        tmp_path: Path,
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Test write_genome_json formats and generates the correct JSON output."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(FeatureLocation(10, 20), type="gene", qualifiers={"gene": ["GlyrA"]})
        record1.features.append(gene_feature1)

        record2 = SeqRecord(Seq("ATGC"), id="record2")
        gene_feature1 = SeqFeature(FeatureLocation(10, 20), type="gene", qualifiers={"gene": ["GlyrB"]})
        record1.features.append(gene_feature1)

        formatted_files_generator.seq_records = [record1, record2]
        genome_path = tmp_path / "genome.json"
        formatted_files_generator.files["genome"] = genome_path
        # pylint: disable=protected-access
        formatted_files_generator._format_genome_data()
        assert genome_path.exists()
        genome_data = get_json(genome_path)
        assert genome_data["species"]["production_name"] == "TEST_prod"
        assert genome_data["added_seq"]["region_name"] == ["record1", "record2"]

    @pytest.mark.dependency(depends=["parse_genbank"])
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._get_codon_table")
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._prepare_location")
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._write_seq_region_json")
    def test_format_seq_region_json(
        self,
        mock_codon_table: Mock,
        mock_org_location: Mock,
        mock_write_seq_json: Mock,
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Check seq_region.json file contains the correct metadata"""
        record = SeqRecord(Seq("ATGC"), id="record", annotations={"topology": "circular"})
        CDS_feature = SeqFeature(
            FeatureLocation(10, 20), type="CDS", qualifiers={"gene": ["GlyrA"], "transl_table": "2"}
        )
        record.features.append(CDS_feature)
        record.annotations["organelle"] = "mitochondrion"
        formatted_files_generator.seq_records = [record]
        # pylint: disable=protected-access
        formatted_files_generator._format_write_seq_json()

        # Call the method to be tested
        mock_codon_table.assert_called_once()
        mock_org_location.assert_called_once()
        mock_write_seq_json.assert_called_once()

    @pytest.mark.dependency(name="format_gff", depends=["parse_genbank"])
    @pytest.mark.parametrize(
        "all_ids, peptides",
        [(["ID1", "ID2", "ID3"], ["pep1", "pep2"])],
    )
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._parse_record")
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._write_genes_gff")
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._write_pep_fasta")
    def test_format_write_genes_gff(
        self,
        mock_write_pep: Mock,
        mock_write_genes: Mock,
        mock_parse_record: Mock,
        all_ids: List[str],
        peptides: List[str],
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Check gene features in GFF3 format are generated as expected."""
        record = SeqRecord(Seq("ATGC"), id="record")
        gene_feature = SeqFeature(FeatureLocation(10, 20), type="gene", qualifiers={"gene": ["GlyrA"]})
        record.features.append(gene_feature)

        formatted_files_generator.seq_records = [record]
        mock_new_record = record

        mock_parse_record.return_value = (
            mock_new_record,  # Mock the new record
            all_ids,  # Mock the list of all IDs
            peptides,  # Mock the list of peptides
        )

        # Call the method to be tested
        # pylint: disable=protected-access
        formatted_files_generator._format_write_genes_gff()

        # Assert that _parse_record was called
        mock_parse_record.assert_called_once()
        mock_write_genes.assert_called_once()
        mock_write_pep.assert_called_once()

        all_ids.append(all_ids[0])

        mock_parse_record.return_value = (
            mock_new_record,  # Mock the new record
            all_ids,  # Mock the list of all IDs
            peptides,  # Mock the list of peptides
        )

        with pytest.raises(GBParseError):
            formatted_files_generator._format_write_genes_gff()

    @pytest.mark.dependency(depends=["parse_genbank", "format_gff"])
    def test_write_genes_gff(
        self,
        tmp_path: Path,
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Test if GFF3 file is generated when there are SeqFeatures present"""

        record = SeqRecord(Seq("ATGC"), id="record")
        gene_feature = SeqFeature(FeatureLocation(10, 20), type="gene", qualifiers={"gene": ["GlyrA"]})
        CDS_feature = SeqFeature(
            FeatureLocation(10, 15), type="CDS", qualifiers={"gene": ["GlyrA"], "transl_table": "2"}
        )
        record.features = [gene_feature, CDS_feature]
        formatted_files_generator.seq_records = [record]

        formatted_files_generator.files["gene_models"] = tmp_path / "genes.gff"
        # pylint: disable=protected-access
        formatted_files_generator._write_genes_gff(formatted_files_generator.seq_records)
        assert (formatted_files_generator.files["gene_models"]).exists()
        for rec in GFF.parse(formatted_files_generator.files["gene_models"]):
            assert rec.id == record.id
            assert len(rec.features) == len(record.features)

    @pytest.mark.dependency(depends=["parse_genbank", "format_gff"])
    def test_write_pep_fasta(
        self,
        tmp_path: Path,
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Test if peptides FASTA file is generated when peptides are identified"""
        record = SeqRecord(Seq("MFLRTQARFFHATTKKM"), id="cds-record")
        CDS_feature = SeqFeature(
            FeatureLocation(10, 20), type="CDS", qualifiers={"gene": ["GlyrA"], "transl_table": "2"}
        )
        record.features.append(CDS_feature)
        formatted_files_generator.files["fasta_pep"] = tmp_path / "pep.fasta"
        # pylint: disable=protected-access
        formatted_files_generator._write_pep_fasta([record])
        assert (tmp_path / "pep.fasta").exists()

        fasta_pep = SeqIO.read((tmp_path / "pep.fasta"), "fasta")
        assert fasta_pep.id == record.id
        assert fasta_pep.seq == record.seq
