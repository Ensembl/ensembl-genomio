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
"""Unit testing of `ensembl.io.genomio.genome_metadata.prepare` module.

Typical usage example::
    $ pytest test_prepare.py

"""

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import Any, Callable, ContextManager, Dict, Optional

from deepdiff import DeepDiff
import pytest

from ensembl.io.genomio.genome_metadata import prepare


# @pytest.mark.dependency(name="test_get_gbff_regions")
@pytest.mark.parametrize(
    "genome_file, gff3_file, output, expectation",
    [
        pytest.param(
            "genbank_genome.json",
            None,
            {
                "assembly": {
                    "provider_name": "GenBank",
                    "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_013436015.2",
                },
            },
            does_not_raise(),
            id="GenBank assembly",
        ),
        pytest.param(
            "refseq_genome.json",
            "fake.gff3",
            {
                "assembly": {
                    "provider_name": "RefSeq",
                    "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004695.1",
                },
                "annotation": {
                    "provider_name": "RefSeq",
                    "provider_url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004695.1",
                },
            },
            does_not_raise(),
            id="RefSeq assembly and annotation",
        ),
        pytest.param(
            "updated_genome.json",
            "fake.gff3",
            {
                "assembly": {"provider_name": "GenBank", "provider_url": None},
                "annotation": {"provider_name": "GenBank", "provider_url": None},
            },
            does_not_raise(),
            id="Provider information already present",
        ),
        pytest.param(
            "cncb_genome.json", None, {}, pytest.raises(prepare.MetadataError), id="Unexpected provider"
        ),
    ],
)
def test_add_provider(
    json_data: Callable[[str], Any],
    genome_file: str,
    gff3_file: Optional[str],
    output: Dict[str, Dict[str, Optional[str]]],
    expectation: ContextManager,
) -> None:
    """Tests the `prepare.add_provider()` method.

    Args:
        json_data: JSON test file parsing fixture.
        genome_file: Genome metadata JSON file.
        gff3_file: GFF3 file.
        output: Expected elements present in the updated genome metadata.
        expectation: Context manager for the expected exception (if any).
    """
    with expectation:
        genome_metadata = json_data(genome_file)
        prepare.add_provider(genome_metadata, gff3_file)
        for section, metadata in output.items():
            for key, value in metadata.items():
                assert genome_metadata[section].get(key, None) == value


@pytest.mark.parametrize(
    "genome_file, output",
    [
        ("genbank_genome.json", 2),
        ("updated_genome.json", 1),
    ],
)
def test_add_assembly_version(json_data: Callable[[str], Any], genome_file: str, output: int) -> None:
    """Tests the `prepare.add_assembly_version()` method.

    Args:
        json_data: JSON test file parsing fixture.
        genome_file: Genome metadata JSON file.
        output: Assembly version expected in the updated genome metadata.
    """
    genome_metadata = json_data(genome_file)
    prepare.add_assembly_version(genome_metadata)
    assert genome_metadata["assembly"]["version"] == output
