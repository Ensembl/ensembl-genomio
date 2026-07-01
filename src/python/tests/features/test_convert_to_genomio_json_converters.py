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
"""Unit testing of convert-to-GenomIO JSON converter registration."""

from pathlib import Path
from unittest.mock import Mock, call, patch

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import converters, repeatmasker


def test_converter_options_defaults_to_no_repeatmasker_consensus_library() -> None:
    """Test converter options default to no RepeatMasker consensus library."""
    assert convert_to_genomio_json.ConverterOptions().repeatmasker_consensus_lib_path is None


def test_repeatmasker_argument_helper_lives_with_repeatmasker_converters() -> None:
    """Test RepeatMasker-specific CLI helpers live in the RepeatMasker module."""
    assert not hasattr(converters, "_add_repeatmasker_common_args")
    assert hasattr(repeatmasker, "_add_repeatmasker_common_args")


def test_converter_registry_contains_supported_logic_names() -> None:
    """Test the explicit converter registry maps supported analysis logic names."""
    assert tuple(convert_to_genomio_json.CONVERTERS_BY_LOGIC_NAME) == (
        "repeatmask_customlib",
        "repeatmask_repbase",
        "trf",
        "repeatdetector",
    )
    assert (
        convert_to_genomio_json.CONVERTERS_BY_LOGIC_NAME["repeatmask_customlib"]
        is convert_to_genomio_json.RepeatMaskerCustomConverter
    )
    assert (
        convert_to_genomio_json.CONVERTERS_BY_LOGIC_NAME["repeatmask_repbase"]
        is convert_to_genomio_json.RepeatMaskerRepbaseConverter
    )
    assert convert_to_genomio_json.CONVERTERS_BY_LOGIC_NAME["trf"] is convert_to_genomio_json.TrfConverter
    assert convert_to_genomio_json.TrfConverter.__module__.endswith(".trf")
    assert convert_to_genomio_json.RepeatMaskerCustomConverter.__module__.endswith(".repeatmasker")
    assert convert_to_genomio_json.RepeatMaskerRepbaseConverter.__module__.endswith(".repeatmasker")


def test_registered_top_level_converters_are_in_cli_order() -> None:
    """Test top-level CLI converters expose the existing subcommand order."""
    # Keep the subject under test on the left for readability.
    assert convert_to_genomio_json.TOP_LEVEL_CONVERTERS == (  # noqa: SIM300
        convert_to_genomio_json.TrfConverter,
        convert_to_genomio_json.RepeatMaskerConverter,
    )


def test_converter_parse_features_uses_tool_specific_parser(
    convert_to_genomio_json_data_dir: Path,
) -> None:
    """Test converter classes call the matching tool parser."""
    input_path = convert_to_genomio_json_data_dir / "trf" / "success_plain_no_params.dat"

    features, consensuses_by_key = convert_to_genomio_json.TrfConverter.parse_features(
        input_path,
        convert_to_genomio_json.ConverterOptions(),
    )

    assert features
    assert consensuses_by_key


@patch("ensembl.io.genomio.features.convert_to_genomio_json.repeatmasker.parse_repeatmasker_output")
def test_repeatmasker_converter_parse_features_delegates_to_parser(
    mock_parse_repeatmasker_output: Mock,
) -> None:
    """Test RepeatMasker converters delegate to the RepeatMasker parser."""
    input_path = Path("input.out")
    consensus_lib_path = Path("consensus.fa")
    options = convert_to_genomio_json.ConverterOptions(repeatmasker_consensus_lib_path=consensus_lib_path)
    expected: converters.ParseFeaturesResult = ([{"seq_region": "chr1"}], {})
    mock_parse_repeatmasker_output.return_value = expected

    assert convert_to_genomio_json.RepeatMaskerCustomConverter.parse_features(input_path, options) == expected
    assert (
        convert_to_genomio_json.RepeatMaskerRepbaseConverter.parse_features(input_path, options) == expected
    )
    assert mock_parse_repeatmasker_output.call_args_list == [
        call(input_path, consensus_lib_path),
        call(input_path, consensus_lib_path),
    ]
