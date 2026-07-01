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
from ensembl.io.genomio.features.convert_to_genomio_json.registry import TOP_LEVEL_CONVERTERS
from ensembl.utils.argparse import ArgumentParser

__all__ = ["parse_args"]


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
    for converter in TOP_LEVEL_CONVERTERS:
        converter.add_parser(subparsers)

    return parser.parse_args(argv)
