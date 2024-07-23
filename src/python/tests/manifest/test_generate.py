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
"""Unit testing of `ensembl.io.genomio.manifest.generate` module."""

from pathlib import Path

import pytest
from pytest import param

from ensembl.io.genomio.manifest.generate import ManifestMaker


_CONTENT_MD5 = "45685e95985e20822fb2538a522a5ccf"
@pytest.mark.parametrize(
    "file_name, expected_name",
    [
        param("", "", id="No files"),
        param("foobar.txt", ""),

        param("gene_models.gff3", "gff3", id="gff3 from name and suffix"),
        param("foobar.gff3", "gff3", id="gff3 from suffix"),

        param("fasta_dna.fasta", "fasta_dna"),
        param("fasta_pep.fasta", "fasta_pep"),
        param("foobar_fasta_dna.fasta", "fasta_dna"),
        param("foobar_fasta_pep.fasta", "fasta_pep"),
        param("foobar.fasta", "", id="Not recognized fasta"),

        param("functional_annotation.json", "functional_annotation"),
        param("foobar_functional_annotation.json", "functional_annotation"),
        param("genome.json", "genome"),
        param("foobar_genome.json", "genome"),
        param("seq_attrib.json", "seq_attrib"),
        param("foobar_seq_attrib.json", "seq_attrib"),
        param("seq_region.json", "seq_region"),
        param("foobar_seq_region.json", "seq_region"),
        param("foobar.json", "", id="Not recognized json"),

        param("foobar.agp", "agp", id="agp from suffix"),
    ],
)
def test_get_files_checksum(tmp_path: Path, file_name: str, expected_name: str) -> None:
    """Tests the `ManifestMaker.get_files_checksum()` method.

    Args:
        tmp_path: Test tmp dir.
        expected_name: Manifest key expected.
        file_name: File to create for the test.

    """
    expected_content = {}
    if file_name:
        with Path(tmp_path / file_name).open("w") as fh:
            fh.write("CONTENT")
        if expected_name:
            expected_content = {expected_name: {"file": file_name, "md5sum": _CONTENT_MD5}}
    maker = ManifestMaker(tmp_path)
    test_files = maker.get_files_checksums()
    assert test_files == expected_content
