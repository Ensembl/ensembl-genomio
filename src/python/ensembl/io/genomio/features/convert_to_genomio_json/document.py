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
"""Create GenomIO JSON documents from parsed repeat feature records."""

import json
from pathlib import Path

from ensembl.io.genomio.features.convert_to_genomio_json.base import file_last_modified_time
from ensembl.io.genomio.features.convert_to_genomio_json.converters import ConverterOptions
from ensembl.io.genomio.features.convert_to_genomio_json.registry import CONVERTERS_BY_LOGIC_NAME

__all__ = ["create_genomio_json"]


def create_genomio_json(  # noqa: PLR0913
    input_path: Path,
    output_path: Path,
    *,
    analysis_logic_name: str,
    analysis_display_label: str,
    analysis_description: str,
    program: str,
    program_version: str,
    source_provider: str,
    is_primary: bool,
    repeatmasker_consensus_lib_path: Path | None = None,
    program_parameters: str | None = None,
) -> None:
    """Create a GenomIO JSON document from feature identification tool output.

    Args:
        input_path: Path to the input file containing repeat masking results (e.g. RepeatMasker .out file).
        output_path: Path to the output JSON file to be created.
        analysis_logic_name: Logic name for the analysis (e.g. "repeatmask_customlib").
        analysis_display_label: Display label for the analysis.
        analysis_description: Description of the analysis (HTML allowed).
        program: Name of the program used to identify the features.
        program_version: Version of the program used to identify the features.
        source_provider: Name of the source provider for the features (e.g. "Ensembl").
        is_primary: Whether these features are primary or secondary annotations.
        repeatmasker_consensus_lib_path: Optional path to a FASTA file containing consensus sequences for the
            RepeatMasker library used.
        program_parameters: Optional parameters supplied to the program used for feature identification.

    Raises:
        ValueError: If an unsupported analysis logic name is provided.

    """
    try:
        converter = CONVERTERS_BY_LOGIC_NAME[analysis_logic_name]
    except KeyError:
        raise ValueError(f"Unsupported analysis logic name: {analysis_logic_name}") from None
    converter_options = ConverterOptions(repeatmasker_consensus_lib_path=repeatmasker_consensus_lib_path)
    features, consensuses_by_key = converter.parse_features(input_path, converter_options)

    analysis: dict[str, str] = {
        "run_date": file_last_modified_time(input_path),
        "logic_name": analysis_logic_name,
        "display_label": analysis_display_label,
        "description": analysis_description,
        "program": program,
        "program_version": program_version,
    }
    if program_parameters is not None:
        analysis["program_parameters"] = program_parameters

    json_doc: dict[str, object] = {
        "analysis": analysis,
        "source": {
            "source_provider": source_provider,
            "is_primary": bool(is_primary),
        },
        "repeat_features": features,
    }

    if consensuses_by_key:
        repeat_consensuses: list[dict[str, str]] = []
        for consensus_key, consensus in consensuses_by_key.items():
            repeat_consensuses.append(
                {
                    "repeat_consensus_key": consensus_key,
                    "repeat_name": consensus.name,
                    "repeat_class": consensus.repeat_class,
                    "repeat_type": consensus.repeat_type,
                    "repeat_consensus": consensus.seq,
                }
            )
        json_doc["repeat_consensus"] = repeat_consensuses

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(json_doc, indent=2) + "\n", encoding="utf-8")
