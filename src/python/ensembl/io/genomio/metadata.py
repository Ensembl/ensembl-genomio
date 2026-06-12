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
"""Create a YAML file with annotation provider acronyms for FTP documentation."""

import argparse
from collections import defaultdict
from importlib.resources import as_file, files
import logging
from pathlib import Path

import yaml

import ensembl.io.genomio.data
from ensembl.io.genomio.utils import get_json
from ensembl.utils import StrPath
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def create_provider_ftp_yaml(output_path: StrPath) -> None:
    """Create a YAML file with provider acronyms for FTP.

    Args:
        output_path: Path to the output YAML file to create.

    """
    source = files(ensembl.io.genomio.data).joinpath("provider_acronyms.json")
    with as_file(source) as json_file:
        provider_acronyms = get_json(json_file)
    # Invert the provider_acronyms mapping to group providers by acronym, as multiple providers may share
    # the same acronym
    yaml_content = defaultdict(list)
    for key, value in provider_acronyms.items():
        yaml_content[value].append(key)
    with Path(output_path).open("w") as file:
        yaml.dump(dict(yaml_content), file)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the FASTA splitting CLI.

    Args:
        argv: Optional argument vector (excluding the program name). If `None`, arguments are read from
            ``sys.argv`` by argparse.

    Returns:
        Parsed argparse namespace with validated options.

    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_dst_path(
        "--output",
        metavar="YAML",
        required=True,
        help="Output YAML file with provider acronyms.",
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Entry point for the FASTA splitting CLI."""
    args = parse_args(argv)
    init_logging_with_args(args)
    try:
        create_provider_ftp_yaml(output_path=args.output)
    except Exception:
        logging.exception(f"Error creating YAML file with provider acronyms for FTP at {args.output}")
        raise
