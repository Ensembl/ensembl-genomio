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
"""Unit testing of `ensembl.io.genomio.fasta.process` module."""
# pylint: disable=too-many-positional-arguments

from contextlib import nullcontext as does_not_raise
import filecmp
from pathlib import Path
from typing import ContextManager, Set

import pytest

import ensembl.io.genomio.fasta.process as FastaProcessing


@pytest.mark.parametrize(
    "input_fasta, input_gbff, pep_mode, expected_output_fasta",
    [
        (
            "input.protein.fa.gz",
            "input.gbff.gz",
            True,
            "output.protein.fa",
        ),
        (
            "input.protein.fa.gz",
            None,
            True,
            "output.protein.fa",
        ),
        (
            "input.protein.fa.gz",
            None,
            False,
            "output.protein.fa",
        ),
        ("input.fna.gz", None, False, "output.nuc.fna"),
        ("input.fna.gz", None, True, "output.nuc.fna"),
    ],
)
def test_fasta_prep(
    tmp_path: Path,
    data_dir: Path,
    input_fasta: str,
    input_gbff: str,
    pep_mode: bool,
    expected_output_fasta: str,
) -> None:
    """Tests the `process.prep_fasta_data()` function.

    Args:
        tmp_path: Where temporary files will be created.
        input_fasta: Name of the fasta file with example input, in the test folder.
        input_gbff: Name of the input GBFF example input, in the test folder.
        pep_mode: Boolean flag to set processing in peptide mode.
        expected_output_fasta: Name of the output fasta file with expected output, in the test folder.
    """

    fasta_input_path = data_dir / input_fasta
    if input_gbff is not None:
        gbff_input_path = data_dir / input_gbff
    else:
        gbff_input_path = input_gbff
    fasta_output_path = tmp_path / f"{input_fasta}.test.fasta"

    FastaProcessing.prep_fasta_data(
        fasta_infile=fasta_input_path,
        genbank_infile=gbff_input_path,
        fasta_outfile=fasta_output_path,
        peptide_mode=pep_mode,
    )

    expected_path = data_dir / expected_output_fasta

    assert filecmp.cmp(fasta_output_path, expected_path)


@pytest.mark.parametrize(
    "input_gbff, excluded_seq_regions, output, expectation",
    [
        (
            "input.gbff.gz",
            set(["LR605957.1"]),
            set(["VWP78966.1", "VWP78967.1", "VWP78968.1"]),
            does_not_raise(),
        ),
        (
            "input.mod.gbff.gz",
            set(["LR605957.1"]),
            set(["VWP78966.1", "VWP78967.1", "VWP78968.1"]),
            pytest.raises(FastaProcessing.FastaParserError),
        ),
    ],
)
def test_exclude_seq_regions(
    data_dir: Path,
    input_gbff: str,
    excluded_seq_regions: Set[str],
    output: Set[str],
    expectation: ContextManager,
) -> None:
    """Tests the `process.get_peptides_to_exclude()` function.

    Args:
        input_gbff: Name of the input GBFF example input, in the test folder.
        excluded_seq_regions: Set of sequence regions to be excluded.
        output: Set of proteins expected to be excluded given excluded_seq_region
        expectation: Context manager for the expected exception, i.e. the test will only pass if that
            exception is raised. Use `~contextlib.nullcontext` if no exception is expected.

    """
    if input_gbff is not None:
        gbff_input_path = data_dir / input_gbff
    else:
        gbff_input_path = input_gbff
    with expectation:
        excluded_proteins = FastaProcessing.get_peptides_to_exclude(gbff_input_path, excluded_seq_regions)
        assert set(excluded_proteins) == set(output)
