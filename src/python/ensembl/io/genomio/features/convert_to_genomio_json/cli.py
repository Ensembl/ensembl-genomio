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
"""Command-line entry point for repeat feature JSON conversion."""

import logging

from ensembl.io.genomio.features.convert_to_genomio_json.args import parse_args
from ensembl.io.genomio.features.convert_to_genomio_json.document import create_genomio_json

__all__ = ["main"]


def main(argv: list[str] | None = None) -> None:
    """Run the JSON conversion command-line entry point.

    Args:
        argv: Optional list of command-line arguments. If `None`, arguments are taken from ``sys.argv``.

    """
    args = parse_args(argv)
    try:
        create_genomio_json(
            input_path=args.input,
            output_path=args.output,
            analysis_logic_name=args.analysis_logic_name,
            analysis_display_label=args.analysis_display_label,
            analysis_description=args.analysis_description,
            program=args.program,
            program_version=args.program_version,
            source_provider=args.source_provider,
            is_primary=args.is_primary,
            repeatmasker_consensus_lib_path=getattr(args, "consensus_lib", None),
            program_parameters=getattr(args, "program_parameters", None),
        )
    except Exception:
        logging.exception(f"Error processing file {args.input}")
        raise
