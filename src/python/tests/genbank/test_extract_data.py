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

    @pytest.mark.parametrize("expected",
                        [
                            ({"NC_000884.1": "mitochondrion"}),
                        ],
                        )
    def test_get_organella(data_dir, expected, formatted_files_generator):
        """Test that organellas are correctly identified."""
        gb_file_path = data_dir/formatted_files_generator.gb_file
        organella = formatted_files_generator._get_organella(gb_file_path)
        assert organella == expected

    def test_parse_genbank(data_dir, formatted_files_generator):
        gb_file_path = data_dir / f"input_file.gb"
        formatted_files_generator.parse_genbank(gb_file_path)
        assert len(formatted_files_generator.seq_records) >= 1

    @pytest.mark.parametrize("type_feature, gene_name, test_qualifiers, expected_id", 
                             [("gene", "AGR90MT_t01", "pseudogene", "TESTAGR90MT_t01" ), 
                              ("tRNA", "AGR90MT_t01", "locus_tag" , "TESTAGR90MT_t01_t1"),
                              ("CDS", "AGR90MT_t01", "translation", "TESTAGR90MT_t01_p1"),
                              ("CDS", "AGR_01", "pseudo", "TESTAGR_01_p1")])
    def test_parse_gene_feat(expected_id: str, gene_name: str, type_feature: str, test_qualifiers: str, formatted_files_generator):
        
        seq_feature = SeqFeature(SimpleLocation(5, 10),type=type_feature, id= gene_name)
        seq_feature.qualifiers[test_qualifiers] = "test_qual"
        # Check the returned feature is as expected
        result_seq_feature, result_seq_id, result_peptide = formatted_files_generator._parse_gene_feat(seq_feature, gene_name)
        if type_feature != "CDS":
            seq_dict = {expected_id: seq_feature}
            assert  result_seq_feature == seq_dict
            assert  result_seq_id == [expected_id]
        else:
            gene_id = formatted_files_generator.prefix + gene_name
            tr_id = gene_id + "_t1"
            assert len(result_seq_feature)> 1
            assert result_seq_id == [tr_id, expected_id]
            if "pseudo" in seq_feature.qualifiers:
                assert result_seq_feature[expected_id].type  == "exon"
            # Peptides aren't always present in the genbank file so we can't guarantee they will exist
            if  "translation" in seq_feature.qualifiers:
                assert len(result_peptide) > 0


    @pytest.mark.parametrize("rna_name, expected_rna_id, expected_gene_id", 
                             [("AGR90MT t01", "AGR90MT t01_t1", "AGR90MT t01"),
                              ("AGR90MT t01 t2", "AGR90MT_t1", "AGR90MT"),
                              ("AGR90MT_t01_t2", "AGR90MT_t01_t2_t1", "AGR90MT_t01_t2")])
    def test_parse_rna_feat(rna_name, expected_gene_id, expected_rna_id, formatted_files_generator):
        seq_feature = SeqFeature(SimpleLocation(5, 10),type="tRNA")
        seq_feature.qualifiers["product"] = [rna_name]
        rna_feature, all_expected_id = formatted_files_generator._parse_rna_feat(seq_feature)
        assert len(rna_feature) == 2
        assert all_expected_id == [formatted_files_generator.prefix+expected_gene_id,formatted_files_generator.prefix+expected_rna_id ]
    
    @pytest.mark.parametrize("gene_id, all_ids, expected_id",  
                             [("gene_name",["gene_name"],"gene_name_2"),
                              ("gene_test",[""],"gene_test")])
    def test_uniquify_id(gene_id, all_ids, expected_id, formatted_files_generator):
        """Test that _uniquify_id adds a version number to an existing ID"""
        new_id = formatted_files_generator._uniquify_id(gene_id,all_ids)
        assert new_id == expected_id

    @pytest.mark.parametrize("organelle, expected_location", [("mitochondrion", "mitochondrial_chromosome")])
    def test_prepare_location_with_supported_organelle( , organelle, formatted_files_generator, expected_location):
        result = formatted_files_generator._prepare_location(organelle)
        assert result == expected_location
    
    @pytest.mark.parametrize("organelle, expected_location", [("miton", "mitochondrial_chromosome")])
    def test_prepare_location_with_unsupported_organelle( , formatted_files_generator, organelle, expected_location):  # An organelle not in the dictionary
        with pytest.raises(UnsupportedData) as exc_info:
            formatted_files_generator._prepare_location(organelle)
        assert str(exc_info.value) == f"Unknown organelle: {organelle}"
        
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

    @pytest.mark.parametrize("type_feature, expected_value", [("gene", None), 
                                    ("mRNA", None), ("CDS", 2), ("CDS", 5)])
    def test_get_codon_table(type_feature, expected_value, data_dir, formatted_files_generator):
        rec = SeqRecord(seq="", id="1JOY", name="EnvZ")
        seq_feature = SeqFeature(type=type_feature, qualifiers={"transl_table": [expected_value]})
        rec.features.append(seq_feature)
        codon_table = formatted_files_generator._get_codon_table(rec)
        assert codon_table == expected_value

    def test_write_seq_region_json(tmp_path, data_dir, formatted_files_generator):
        out_dir = tmp_path
        formatted_files_generator.extract_gb(out_dir)
        output= tmp_path/ "seq_region.json"
        expected_path = data_dir/ f"output_seq_region.json"
        assert filecmp.cmp(output, expected_path)

    @pytest.mark.parametrize("organelle, expected_location", [("mitochondrion", "mitochondrial_chromosome")])
    def test_prepare_location_with_supported_organelle( organelle, formatted_files_generator, expected_location):
        result = formatted_files_generator._prepare_location(organelle)
        assert result == expected_location
    
    @pytest.mark.parametrize("organelle, expected_location", [("miton", "mitochondrial_chromosome")])
    def test_prepare_location_with_unsupported_organelle(formatted_files_generator, organelle, expected_location):  # An organelle not in the dictionary
        with pytest.raises(UnsupportedData) as exc_info:
            formatted_files_generator._prepare_location(organelle)
        assert str(exc_info.value) == f"Unknown organelle: {organelle}"



