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
    $ pytest test_extract_data.py

"""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest
import json

from pathlib import Path

from ensembl.io.genomio.genbank.extract_data import FormattedFilesGenerator, GBParseError, UnsupportedData


class TestFormattedFilesGenerator:
    @pytest.fixture
    def formatted_files_generator(self, data_dir):
        """Call the function `FormattedFilesGenerator` with set parameters"""
        prod_name = "TEST_prod"
        gb_file = "input_file.gb"
        gb_file_path = data_dir / gb_file
        prefix = "TEST"
        return FormattedFilesGenerator(prod_name, gb_file_path, prefix)

    @pytest.mark.parametrize(
        "type_feature, gene_name, test_qualifiers, expected_id",
        [
            ("gene", "AGR90MT_t01", "pseudogene", "TESTAGR90MT_t01"),
            ("tRNA", "AGR90MT_t01", "locus_tag", "TESTAGR90MT_t01_t1"),
            ("CDS", "AGR90MT_t01", "translation", "TESTAGR90MT_t01_p1"),
            ("CDS", "AGR_01", "pseudo", "TESTAGR_01_p1"),
        ],
    )
    def test_parse_gene_feat(
        self,
        expected_id: str,
        gene_name: str,
        type_feature: str,
        test_qualifiers: str,
        formatted_files_generator,
    ) -> None:
        """Test for a successful parsing of `_parse_gene_feat()` method"""
        seq_feature = SeqFeature(SimpleLocation(5, 10), type=type_feature, id=gene_name)
        seq_feature.qualifiers[test_qualifiers] = "test_qual"
        # Check the returned feature is as expected
        result_seq_feature, result_seq_id, result_peptide = formatted_files_generator._parse_gene_feat(
            seq_feature, gene_name
        )

        gene_id = formatted_files_generator.prefix + gene_name
        if type_feature == "gene":
            seq_dict = {expected_id: seq_feature}
            assert result_seq_feature == seq_dict
            assert result_seq_id == [expected_id]
            # Verify the transformation results
            assert "gene" not in seq_feature.qualifiers
            assert "locus_tag" not in seq_feature.qualifiers
            assert "Name" in seq_feature.qualifiers

        if type_feature in ("tRNA", "rRNA"):
            assert "gene" not in seq_feature.qualifiers
            assert "locus_tag" not in seq_feature.qualifiers
            assert seq_feature.qualifiers["Parent"] == gene_id
            assert result_seq_id == [expected_id]
        
        if type_feature == "CDS":
            tr_id = gene_id + "_t1"
            assert len(result_seq_feature) > 1
            assert seq_feature.qualifiers["Parent"] == tr_id
            assert result_seq_id == [tr_id, expected_id]
            if "pseudo" in seq_feature.qualifiers:
                assert result_seq_feature[expected_id].type == "exon"
            # Peptides aren't always present in the genbank file so we can't guarantee they will exist
            if "translation" in seq_feature.qualifiers:
                assert len(result_peptide) > 0           

    @pytest.mark.parametrize(
        "rna_name, expected_rna_id, expected_gene_id",
        [
            ("AGR90MT t01", "AGR90MT t01_t1", "AGR90MT t01"),
            ("AGR90MT t01 t2", "AGR90MT_t1", "AGR90MT"),
            ("AGR90MT_t01_t2", "AGR90MT_t01_t2_t1", "AGR90MT_t01_t2"),
        ],
    )
    def test_parse_rna_feat(
        self, expected_gene_id: str, expected_rna_id: str, formatted_files_generator, rna_name: str
    ) -> None:
        """Test for a successful parsing of `_parse_rna_feat()` method"""
        seq_feature = SeqFeature(SimpleLocation(5, 10), type="tRNA")
        seq_feature.qualifiers["product"] = [rna_name]
        rna_feature, all_expected_id = formatted_files_generator._parse_rna_feat(seq_feature)
        assert len(rna_feature) == 2
        assert all_expected_id == [
            formatted_files_generator.prefix + expected_gene_id,
            formatted_files_generator.prefix + expected_rna_id,
        ]

    @pytest.mark.parametrize(
        "gene_id, all_ids, expected_id",
        [("gene_name", ["gene_name"], "gene_name_2"), ("gene_test", [""], "gene_test")],
    )
    def test_uniquify_id(
        self, all_ids: str, expected_id: str, formatted_files_generator, gene_id: str
    ) -> None:
        """Test that _uniquify_id adds a version number to an existing ID"""
        new_id = formatted_files_generator._uniquify_id(gene_id, all_ids)
        assert new_id == expected_id

    @pytest.mark.parametrize("organelle, expected_location", [("mitochondrion", "mitochondrial_chromosome")])
    def test_prepare_location_with_supported_organelle(
        self, expected_location: str, formatted_files_generator, organelle: str
    ) -> None:
        """Test that organelle location is present in the allowed types"""
        result = formatted_files_generator._prepare_location(organelle)
        assert result == expected_location

    @pytest.mark.parametrize("organelle", [("miton")])
    def test_prepare_location_with_unsupported_organelle(
        self, formatted_files_generator, organelle: str
    ) -> None:
        """Test that organelle location if not identifies throws an error"""
        # An organelle not in the dictionary
        with pytest.raises(UnsupportedData) as exc_info:
            formatted_files_generator._prepare_location(organelle)
        assert str(exc_info.value) == f"Unknown organelle: {organelle}"

    @pytest.mark.parametrize(
        "type_feature, expected_value", [("gene", None), ("mRNA", None), ("CDS", 2), ("CDS", 5)]
    )
    def test_get_codon_table(self, expected_value: str, formatted_files_generator, type_feature: str) -> None:
        """Test that `get_number_of_codons` returns correct value based on feature type and qualifier"""
        rec = SeqRecord(seq="", id="1JOY", name="EnvZ")
        seq_feature = SeqFeature(type=type_feature, qualifiers={"transl_table": [expected_value]})
        rec.features.append(seq_feature)
        codon_table = formatted_files_generator._get_codon_table(rec)
        assert codon_table == expected_value

    def test_parse_record_with_unsupported_feature(self, formatted_files_generator):
        record = SeqRecord(Seq("ATGC"))
        unsupported_feature = SeqFeature(SimpleLocation(5, 10), type="gene")
        record.features.append(unsupported_feature)

        with pytest.raises(GBParseError):
            formatted_files_generator._parse_record(record)

    def test_format_genome_json_with_production_name(self, formatted_files_generator, tmp_path, caplog):
        """Test write_genome_json generates correct json output."""
        record1 = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrA"]})        
        record1.features.append(gene_feature1)

        record2 = SeqRecord(Seq("ATGC"), id="record2")
        gene_feature1 = SeqFeature(SimpleLocation(10,20), type="gene", qualifiers= {"gene": ["GlyrB"]})        
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

