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
    $ pytest test_extract_data_seq.py

"""
from unittest.mock import Mock, patch
from typing import List
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pytest


from ensembl.io.genomio.genbank.extract_data import FormattedFilesGenerator, GBParseError, UnsupportedData


class TestFormattedFilesGenerator:
    """Test if all the internal methods of `FormattedFilesGenerator` are giving the correct output"""

    prod_name = "TEST_prod"
    gb_file = "input_file.gb"
    prefix = "TEST"

    @pytest.fixture
    def formatted_files_generator(self, data_dir):
        """Call the function `FormattedFilesGenerator` with set parameters"""
        gb_file = self.gb_file
        gb_file_path = data_dir / gb_file
        return FormattedFilesGenerator(self.prod_name, gb_file_path, self.prefix)

    @pytest.mark.dependency(name="parse_record")
    @pytest.mark.parametrize(
        "type_feature, gene_name, expected_name, expected_id",
        [
            ("rRNA", "locus_tag", "gene1", "TESTGlyrA"),
            ("tRNA", "gene", "gene2", "TESTGlyrA"),
            ("mRNA", "gene", "gene3", "TESTGlyrA"),
        ],
    )
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._parse_gene_feat")
    @patch("ensembl.io.genomio.genbank.extract_data.FormattedFilesGenerator._parse_rna_feat")
    def test_parse_record(
        self,
        mock_parse_rna_feat: Mock,
        mock_parse_gene_feat: Mock,
        expected_name: str,
        expected_id: str,
        formatted_files_generator: FormattedFilesGenerator,
        gene_name: str,
        type_feature: str,
    ):
        """Test to parse the features correctly"""
        record = SeqRecord(Seq("ATGC"), id="record1")
        gene_feature = SeqFeature(FeatureLocation(10, 20), type="gene", qualifiers={gene_name: expected_name})
        rna_feature = SeqFeature(FeatureLocation(10, 15), type=type_feature)
        CDS_feature = SeqFeature(
            FeatureLocation(10, 20), type="CDS", qualifiers={gene_name: "GlyrA", "transl_table": "2"}
        )
        record.features = [gene_feature, rna_feature, CDS_feature]
        mock_peptides: List = []

        assert gene_feature.qualifiers[gene_name] == expected_name
        gene_feature_feat = {expected_id: gene_feature}
        mock_parse_gene_feat.return_value = (
            gene_feature_feat,  # Mock the new record
            gene_name,  # Mock the list of all IDs
            mock_peptides,  # Mock the list of peptides
        )
        rna_feature_feat = {expected_id: rna_feature}
        mock_parse_rna_feat.return_value = (
            rna_feature_feat,  # Mock the new record
            gene_name,  # Mock the list of peptides
        )
        # pylint: disable=protected-access
        formatted_files_generator._parse_record(record)
        if gene_feature.qualifiers[gene_name]:
            mock_parse_gene_feat.assert_called()

        if rna_feature.type in ("tRNA", "rRNA"):
            mock_parse_rna_feat.assert_called()

    def test_parse_record_with_unsupported_feature(self, formatted_files_generator: FormattedFilesGenerator):
        """Raises error then there is unspported feature"""
        record = SeqRecord(Seq("ATGC"))
        unsupported_feature = SeqFeature(FeatureLocation(5, 10), type="gene")
        record.features.append(unsupported_feature)

        with pytest.raises(GBParseError):
            # pylint: disable=protected-access
            formatted_files_generator._parse_record(record)

    @pytest.mark.dependency(depends=["parse_record"])
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
        formatted_files_generator: FormattedFilesGenerator,
    ) -> None:
        """Test for a successful parsing of features in `_parse_gene_feat()` method"""
        seq_feature = SeqFeature(FeatureLocation(5, 10), type=type_feature, id=gene_name)
        seq_feature.qualifiers[test_qualifiers] = "test_qual"
        # Check the returned feature is as expected
        # pylint: disable=protected-access
        result_seq_feature, result_seq_id, result_peptide = formatted_files_generator._parse_gene_feat(
            seq_feature, gene_name
        )

        gene_id = self.prefix + gene_name
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
    @pytest.mark.dependency(name="rna_parse", depends=["parse_record"])
    def test_parse_rna_feat(
        self,
        expected_gene_id: str,
        expected_rna_id: str,
        formatted_files_generator: FormattedFilesGenerator,
        rna_name: str,
    ) -> None:
        """Test for a successful parsing of transcript features `_parse_rna_feat()` method"""
        seq_feature = SeqFeature(FeatureLocation(5, 10), type="tRNA")
        seq_feature.qualifiers["product"] = [rna_name]
        # pylint: disable=protected-access
        rna_feature, all_expected_id = formatted_files_generator._parse_rna_feat(seq_feature)
        assert len(rna_feature) == 2
        assert all_expected_id == [
            self.prefix + expected_gene_id,
            self.prefix + expected_rna_id,
        ]

    @pytest.mark.dependency(depends=["rna_parse"])
    @pytest.mark.parametrize(
        "gene_id, all_ids, expected_id",
        [("gene_name", ["gene_name"], "gene_name_2"), ("gene_test", [""], "gene_test")],
    )
    def test_uniquify_id(
        self,
        all_ids: List[str],
        expected_id: str,
        formatted_files_generator: FormattedFilesGenerator,
        gene_id: str,
    ) -> None:
        """Test that _uniquify_id adds a version number to an existing ID"""
        # pylint: disable=protected-access
        new_id = formatted_files_generator._uniquify_id(gene_id, all_ids)
        assert new_id == expected_id

    @pytest.mark.parametrize("organelle, expected_location", [("mitochondrion", "mitochondrial_chromosome")])
    def test_prepare_location_with_supported_organelle(
        self, expected_location: str, formatted_files_generator, organelle: str
    ) -> None:
        """Test that organelle location is present in the allowed types"""
        # pylint: disable=protected-access
        result = formatted_files_generator._prepare_location(organelle)
        assert result == expected_location

    @pytest.mark.parametrize("organelle", [("miton")])
    def test_prepare_location_with_unsupported_organelle(
        self, formatted_files_generator: FormattedFilesGenerator, organelle: str
    ) -> None:
        """Test that organelle location if not identifies throws an error"""
        # An organelle not in the dictionary
        with pytest.raises(UnsupportedData) as exc_info:
            # pylint: disable=protected-access
            formatted_files_generator._prepare_location(organelle)
        assert str(exc_info.value) == f"Unknown organelle: {organelle}"

    @pytest.mark.parametrize(
        "type_feature, expected_value", [("gene", None), ("mRNA", None), ("CDS", 2), ("CDS", 5)]
    )
    def test_get_codon_table(
        self, expected_value: str, formatted_files_generator: FormattedFilesGenerator, type_feature: str
    ) -> None:
        """Test that `get_number_of_codons` returns correct value based on feature type and qualifier"""
        rec = SeqRecord(seq="", id="1JOY", name="EnvZ")
        seq_feature = SeqFeature(type=type_feature, qualifiers={"transl_table": [expected_value]})
        rec.features.append(seq_feature)
        # pylint: disable=protected-access
        codon_table = formatted_files_generator._get_codon_table(rec)
        assert codon_table == expected_value
