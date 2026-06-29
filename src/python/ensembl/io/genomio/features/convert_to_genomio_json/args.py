# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Parse command-line arguments for repeat feature JSON conversion."""

import argparse

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser

__all__ = ["parse_args"]


def _add_common_arguments(subparser: ArgumentParser) -> None:
    """Add arguments shared by all supported analysis subcommands."""
    subparser.add_argument_src_path(
        "--input",
        metavar="IN",
        required=True,
        help="Input file to be converted.",
    )
    subparser.add_argument_dst_path("--output", metavar="JSON", required=True, help="Output JSON path.")
    subparser.add_argument(
        "--program-version",
        required=True,
        help="Version of the program used to identify the features.",
    )
    subparser.add_argument(
        "--program-parameters",
        default=argparse.SUPPRESS,
        help="Parameters supplied to the program used for feature identification.",
    )
    subparser.add_argument(
        "--source-provider",
        default="Ensembl",
        help="Source provider for the features (default: Ensembl).",
    )
    subparser.add_argument(
        "--is-primary",
        action="store_true",
        help="Whether the source provider is the primary source for these features.",
    )
    subparser.add_log_arguments()


def _add_repeatmasker_common_args(subparser: ArgumentParser) -> None:
    """Add arguments shared by RepeatMasker subcommands."""
    _add_common_arguments(subparser)
    subparser.add_argument_src_path(
        "--consensus-lib",
        metavar="RM_LIB",
        default=argparse.SUPPRESS,
        help="FASTA file containing consensus sequences for the RepeatMasker library used.",
    )


def _add_trf_parser(subparsers: argparse._SubParsersAction) -> None:  # noqa: SLF001
    """Add the TRF subcommand parser."""
    trf_parser = subparsers.add_parser("trf", help="Convert TRF output to GenomIO JSON.")
    _add_common_arguments(trf_parser)
    trf_parser.set_defaults(
        analysis_logic_name="trf",
        analysis_display_label="Tandem repeats (TRF)",
        analysis_description=(
            '<a rel="external" href="https://tandem.bu.edu/trf/trf.html">Tandem Repeats Finder</a> '
            "locates adjacent copies of a pattern of nucleotides."
        ),
        program="trf",
        repeatmasker_consensus_lib_path=None,
    )


def _add_repeatmasker_parser(subparsers: argparse._SubParsersAction) -> None:  # noqa: SLF001
    """Add the RepeatMasker subcommand parser."""
    repeatmasker_parser = subparsers.add_parser(
        "repeatmasker", help="Convert RepeatMasker output to GenomIO JSON."
    )
    repeatmasker_parser.set_defaults(program="RepeatMasker")
    repeatmasker_subparsers = repeatmasker_parser.add_subparsers(dest="repeatmasker_mode", required=True)

    custom_parser = repeatmasker_subparsers.add_parser(
        "custom",
        help="Convert RepeatMasker output generated using a custom library.",
    )
    _add_repeatmasker_common_args(custom_parser)
    custom_parser.set_defaults(
        analysis_logic_name="repeatmask_customlib",
        analysis_display_label="Repeats: Custom library",
        analysis_description=(
            'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, '
            "using a custom library of <em>ab initio</em> repeat profiles for this species."
        ),
    )

    repbase_parser = repeatmasker_subparsers.add_parser(
        "repbase",
        help="Convert RepeatMasker output geneterated using Repbase.",
    )
    _add_repeatmasker_common_args(repbase_parser)
    repbase_parser.set_defaults(
        analysis_logic_name="repeatmask_repbase",
        analysis_display_label="Repeats: Repbase",
        analysis_description=(
            'Repeats identified by <a rel="external" href="http://www.repeatmasker.org">RepeatMasker</a>, '
            'using the <a rel="external" href="http://www.girinst.org/repbase/">Repbase</a> library of '
            "repeat profiles."
        ),
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for JSON conversion.

    Args:
        argv: Optional list of command-line arguments. If `None`, arguments are
            taken from ``sys.argv``.

    Returns:
        Parsed command-line arguments.

    """
    parser = ArgumentParser(description="Constructs a GenomIO JSON document from feature output.")
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)

    subparsers = parser.add_subparsers(dest="tool", required=True)
    _add_trf_parser(subparsers)
    _add_repeatmasker_parser(subparsers)

    return parser.parse_args(argv)
