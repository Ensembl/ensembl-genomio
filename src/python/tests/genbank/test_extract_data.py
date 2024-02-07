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
    FormattedFilesGenerator,
)

class TestFormattedFilesGenerator:

    prod_name= "TEST_prod"
    gb_file = "input_file.gb"
    prefix = "TEST"

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

    @pytest.mark.parametrize("type_feature, gene_name, test_qualifiers, expected_id", 
                             [("gene", "AGR90MT_t01", "pseudogene", "TESTAGR90MT_t01" ), 
                              ("tRNA", "AGR90MT_t01", "locus_tag" , "TESTAGR90MT_t01_t1"),
                              ("CDS", "AGR90MT_t01", "translation", "TESTAGR90MT_t01_p1"),
                              ("CDS", "AGR_01", "pseudo", "TESTAGR_01_p1")])
    def test_parse_gene_feat(self, expected_id: str, gene_name: str, type_feature: str, test_qualifiers: str):
        
        gb = FormattedFilesGenerator(self.prod_name,self.gb_file, prefix="TEST")
        
        seq_feature = SeqFeature(SimpleLocation(5, 10),type=type_feature, id= gene_name)
        seq_feature.qualifiers[test_qualifiers] = "test_unkown"
        # Check the returned feature is as expected
        result_seq_feature, result_seq_id, result_peptide = gb._parse_gene_feat(seq_feature, gene_name)
        if type_feature != "CDS":
            seq_dict = {expected_id: seq_feature}
            assert  result_seq_feature == seq_dict
            assert  result_seq_id == [expected_id]
        else:
            gene_id = self.prefix + gene_name
            tr_id = gene_id + "_t1"
            assert len(result_seq_feature)> 1
            assert result_seq_id == [tr_id, expected_id]
            if "pseudo" in seq_feature.qualifiers:
                assert result_seq_feature[expected_id].type  == "exon"
            # Peptides aren't always present in the genbank file so we can't guarantee they will exist
            if  "translation" in seq_feature.qualifiers:
                assert len(result_peptide) > 0
        
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

    @pytest.mark.parametrize("type_feature, expected_value", [("gene", None), 
                                    ("mRNA", None), ("CDS", 2), ("CDS", 5)])
    def test_get_codon_table(self, type_feature, expected_value, data_dir):
        gb_file_path = data_dir / self.gb_file
        gb = FormattedFilesGenerator(self.prod_name,gb_file_path)
        rec = SeqRecord(seq="", id="1JOY", name="EnvZ")
        seq_feature = SeqFeature(type=type_feature, qualifiers={"transl_table": [expected_value]})
        rec.features.append(seq_feature)
        codon_table = gb._get_codon_table(rec)
        assert codon_table == expected_value

    def test_write_seq_region_json(self, tmp_path, data_dir):
        gb_file_path = data_dir / self.gb_file
        out_dir = tmp_path
        gb = FormattedFilesGenerator(self.prod_name,gb_file_path)
        gb.extract_gb(out_dir)
        output= tmp_path/ "seq_region.json"
        expected_path = data_dir/ f"output_seq_region.json"
        assert filecmp.cmp(output, expected_path)

