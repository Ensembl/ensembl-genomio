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
"""Construct a seq_region metadata file from INSDC files."""

from pathlib import Path

import ensembl.io.genomio
from ensembl.io.genomio.utils import get_json, print_json
from ensembl.io.genomio.seq_region.collection import SeqCollection
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args
from ensembl.utils import StrPath


def prepare_seq_region_metadata(
    genome_file: StrPath,
    report_file: StrPath,
    dst_file: StrPath,
    *,
    gbff_file: StrPath | None = None,
    to_exclude: list[str] | None = None,
    mock_run: bool = False,
) -> None:
    """Prepares the sequence region metadata found in the INSDC/RefSeq report and GBFF files.

    The sequence region information is loaded from both sources and combined. Elements are added/excluded
    as requested, and the final sequence region metadata is dumped in a JSON file that follows the schema
    defined in "src/python/ensembl/io/genomio/data/schemas/seq_region.json".

    Args:
        genome_file: Genome metadata JSON file path.
        report_file: INSDC/RefSeq sequences report file path to parse.
        gbff_file: INSDC/RefSeq GBFF file path to parse.
        dst_file: JSON file output for the processed sequence regions JSON.
        to_exclude: Sequence region names to exclude.
        mock_run: Do not call external taxonomy service.

    """
    genome_data = get_json(genome_file)
    dst_file = Path(dst_file)
    is_refseq = genome_data["assembly"]["accession"].startswith("GCF_")

    seqs = SeqCollection(mock=mock_run)
    seqs.from_report(Path(report_file), is_refseq)
    if gbff_file:
        seqs.from_gbff(Path(gbff_file))

    # Exclude seq_regions from a list
    if to_exclude:
        seqs.remove(to_exclude)

    # Add translation and mitochondrial codon tables
    seqs.add_translation_table()
    seqs.add_mitochondrial_codon_table(genome_data["species"]["taxonomy_id"])

    # Print out the file
    print_json(dst_file, seqs.to_list())


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description="Construct a sequence region metadata file from INSDC files.")
    parser.add_argument_src_path("--genome_file", required=True, help="Genome metadata JSON file")
    parser.add_argument_src_path(
        "--report_file", required=True, help="INSDC/RefSeq sequences report file to parse"
    )
    parser.add_argument_src_path("--gbff_file", help="INSDC/RefSeq GBFF file to parse")
    parser.add_argument_dst_path(
        "--dst_file", default="seq_region.json", help="Output JSON file for the processed sequence regions"
    )
    parser.add_argument(
        "--to_exclude", nargs="*", metavar="SEQ_REGION_NAME", help="Sequence region names to exclude"
    )
    parser.add_argument("--mock_run", action="store_true", help="Do not call external APIs")
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    prepare_seq_region_metadata(
        genome_file=args.genome_file,
        report_file=args.report_file,
        dst_file=args.dst_file,
        gbff_file=args.gbff_file,
        to_exclude=args.to_exclude,
        mock_run=args.mock_run,
    )
