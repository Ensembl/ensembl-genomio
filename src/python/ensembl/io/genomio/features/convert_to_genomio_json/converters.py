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
"""Converter classes for feature output tools."""

import argparse
from dataclasses import dataclass
from pathlib import Path

from ensembl.io.genomio.features.convert_to_genomio_json.base import Consensus
from ensembl.utils.argparse import ArgumentParser

__all__ = [
    "ConverterOptions",
    "FeatureConverter",
    "ParseFeaturesResult",
]


ParseFeaturesResult = tuple[list[dict[str, object]], dict[str, Consensus]]


@dataclass(frozen=True)
class ConverterOptions:
    """Tool-specific options supplied to feature converters."""

    repeatmasker_consensus_lib_path: Path | None = None


class FeatureConverter:
    """Base class for feature output converters."""

    analysis_logic_name: str | None = None
    command: str | None = None

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add this converter's CLI parser."""
        raise NotImplementedError

    @classmethod
    def parse_features(
        cls,
        input_path: Path,
        options: ConverterOptions | None = None,
    ) -> ParseFeaturesResult:
        """Parse features and consensus records from a tool output file."""
        raise NotImplementedError


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
