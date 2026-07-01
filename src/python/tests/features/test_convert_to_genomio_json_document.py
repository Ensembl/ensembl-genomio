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
"""Unit testing of GenomIO JSON document assembly."""

import json
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json


@patch(
    "ensembl.io.genomio.features.convert_to_genomio_json.repeatmasker."
    "RepeatMaskerCustomConverter.parse_features"
)
def test_create_genomio_json_uses_repeatmasker_parser_output(
    mock_parse_features: Mock,
    tmp_path: Path,
) -> None:
    """Test JSON creation assembles mocked RepeatMasker parser output.

    Args:
        mock_parse_features: Mock for ``RepeatMaskerCustomConverter.parse_features()``.
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = tmp_path / "input.out"
    input_path.write_text("parser input", encoding="utf-8")
    output_path = tmp_path / "out.json"
    consensus_lib = Path("consensus.fa")
    expected_features = [
        {
            "seq_region": "chr1",
            "seq_region_start": 10,
            "seq_region_end": 20,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 11,
            "repeat_consensus": "consensus-key",
            "score": 100.0,
            "attributes": {"repeatmasker_repeat_type": "Alu"},
        }
    ]
    expected_consensus = convert_to_genomio_json.Consensus(
        name="Foo",
        repeat_class="SINE",
        repeat_type="Type I Transposons/SINE",
        seq="acgt",
    )
    mock_parse_features.return_value = (expected_features, {"consensus-key": expected_consensus})

    convert_to_genomio_json.create_genomio_json(
        input_path=input_path,
        repeatmasker_consensus_lib_path=consensus_lib,
        output_path=output_path,
        analysis_logic_name="repeatmask_customlib",
        analysis_display_label="Repeats: Custom library",
        analysis_description="desc",
        program="RepeatMasker",
        program_version="4.1.5",
        program_parameters="-nolow -gccalc -q",
        source_provider="Ensembl",
        is_primary=False,
    )

    mock_parse_features.assert_called_once_with(
        input_path,
        convert_to_genomio_json.ConverterOptions(repeatmasker_consensus_lib_path=consensus_lib),
    )

    doc = json.loads(output_path.read_text(encoding="utf-8"))

    assert doc["analysis"]["logic_name"] == "repeatmask_customlib"
    assert doc["analysis"]["display_label"] == "Repeats: Custom library"
    assert doc["analysis"]["description"] == "desc"
    assert doc["analysis"]["program"] == "RepeatMasker"
    assert doc["analysis"]["program_version"] == "4.1.5"
    assert doc["analysis"]["program_parameters"] == "-nolow -gccalc -q"
    assert doc["analysis"]["run_date"].endswith("Z")

    assert doc["source"] == {"source_provider": "Ensembl", "is_primary": False}
    assert doc["repeat_features"] == expected_features
    assert doc["repeat_consensus"] == [
        {
            "repeat_consensus_key": "consensus-key",
            "repeat_name": "Foo",
            "repeat_class": "SINE",
            "repeat_type": "Type I Transposons/SINE",
            "repeat_consensus": "acgt",
        }
    ]


@patch("ensembl.io.genomio.features.convert_to_genomio_json.trf.TrfConverter.parse_features")
def test_create_genomio_json_uses_trf_parser_output(
    mock_parse_features: Mock,
    tmp_path: Path,
) -> None:
    """Test JSON creation assembles mocked TRF parser output.

    Args:
        mock_parse_features: Mock for ``TrfConverter.parse_features()``.
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = tmp_path / "input.out"
    input_path.write_text("parser input", encoding="utf-8")
    output_path = tmp_path / "out.json"
    expected_features = [
        {
            "seq_region": "chrA",
            "seq_region_start": 104,
            "seq_region_end": 109,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 2,
            "repeat_consensus": "trf-key",
            "score": 42.0,
            "attributes": {"period_size": 2},
        }
    ]
    expected_consensus = convert_to_genomio_json.Consensus(
        name="trf",
        repeat_class="trf",
        repeat_type="Tandem repeats",
        seq="AT",
    )
    mock_parse_features.return_value = (expected_features, {"trf-key": expected_consensus})

    convert_to_genomio_json.create_genomio_json(
        input_path=input_path,
        repeatmasker_consensus_lib_path=None,
        output_path=output_path,
        analysis_logic_name="trf",
        analysis_display_label="Repeats: Custom library",
        analysis_description="desc",
        program="TRF",
        program_version="4.10.0",
        program_parameters="-h",
        source_provider="Ensembl",
        is_primary=False,
    )

    mock_parse_features.assert_called_once_with(input_path, convert_to_genomio_json.ConverterOptions())

    doc = json.loads(output_path.read_text(encoding="utf-8"))

    assert doc["analysis"]["logic_name"] == "trf"
    assert doc["analysis"]["program"] == "TRF"
    assert doc["analysis"]["program_version"] == "4.10.0"
    assert doc["analysis"]["program_parameters"] == "-h"
    assert doc["source"] == {"source_provider": "Ensembl", "is_primary": False}
    assert doc["repeat_features"] == expected_features
    assert doc["repeat_consensus"] == [
        {
            "repeat_consensus_key": "trf-key",
            "repeat_name": "trf",
            "repeat_class": "trf",
            "repeat_type": "Tandem repeats",
            "repeat_consensus": "AT",
        }
    ]


@patch("ensembl.io.genomio.features.convert_to_genomio_json.red.RedConverter.parse_features")
def test_create_genomio_json_uses_red_parser_output(
    mock_parse_features: Mock,
    tmp_path: Path,
) -> None:
    """Test JSON creation assembles mocked Red parser output.

    Args:
        mock_parse_features: Mock for ``RedConverter.parse_features()``.
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = tmp_path / "input.rpt"
    input_path.write_text("parser input", encoding="utf-8")
    output_path = tmp_path / "out.json"
    expected_features = [
        {
            "seq_region": "chr1",
            "seq_region_start": 10,
            "seq_region_end": 20,
            "seq_region_strand": "+",
            "repeat_start": 1,
            "repeat_end": 11,
            "repeat_consensus": "red-key",
        }
    ]
    expected_consensus = convert_to_genomio_json.Consensus(
        name="repeatdetector",
        repeat_class="repeatdetector",
        repeat_type="repeatdetector",
        seq="N",
    )
    mock_parse_features.return_value = (expected_features, {"red-key": expected_consensus})

    convert_to_genomio_json.create_genomio_json(
        input_path=input_path,
        repeatmasker_consensus_lib_path=None,
        output_path=output_path,
        analysis_logic_name="repeatdetector",
        analysis_display_label="Repeats: Red",
        analysis_description="desc",
        program="Red",
        program_version="2.0",
        program_parameters="-gnm genome.fa",
        source_provider="Custom",
        is_primary=True,
    )

    mock_parse_features.assert_called_once_with(input_path, convert_to_genomio_json.ConverterOptions())

    doc = json.loads(output_path.read_text(encoding="utf-8"))

    assert doc["analysis"]["logic_name"] == "repeatdetector"
    assert doc["analysis"]["program"] == "Red"
    assert doc["analysis"]["program_version"] == "2.0"
    assert doc["analysis"]["program_parameters"] == "-gnm genome.fa"
    assert doc["source"] == {"source_provider": "Custom", "is_primary": True}
    assert doc["repeat_features"] == expected_features
    assert doc["repeat_consensus"] == [
        {
            "repeat_consensus_key": "red-key",
            "repeat_name": "repeatdetector",
            "repeat_class": "repeatdetector",
            "repeat_type": "repeatdetector",
            "repeat_consensus": "N",
        }
    ]


def test_create_genomio_json_rejects_unsupported_logic_name(
    convert_to_genomio_json_data_dir: Path, tmp_path: Path
) -> None:
    """Test JSON creation rejects unsupported analysis logic names.

    Args:
        convert_to_genomio_json_data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.

    """
    with pytest.raises(ValueError, match=r"Unsupported analysis logic name"):
        convert_to_genomio_json.create_genomio_json(
            input_path=convert_to_genomio_json_data_dir / "create_json" / "basic.out",
            repeatmasker_consensus_lib_path=None,
            output_path=tmp_path / "out.json",
            analysis_logic_name="unsupported",
            analysis_display_label="label",
            analysis_description="desc",
            program="RepeatMasker",
            program_version="4.1.5",
            program_parameters=None,
            source_provider="Ensembl",
            is_primary=False,
        )


@patch(
    "ensembl.io.genomio.features.convert_to_genomio_json.repeatmasker."
    "RepeatMaskerCustomConverter.parse_features"
)
def test_create_genomio_json_omits_program_parameters_when_none(
    mock_parse_features: Mock,
    tmp_path: Path,
) -> None:
    """Test that program parameters are omitted when no value is provided.

    Args:
        mock_parse_features: Mock for ``RepeatMaskerCustomConverter.parse_features()``.
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = tmp_path / "input.out"
    input_path.write_text("parser input", encoding="utf-8")
    output_path = tmp_path / "out.json"
    mock_parse_features.return_value = ([], {})

    convert_to_genomio_json.create_genomio_json(
        input_path=input_path,
        repeatmasker_consensus_lib_path=None,
        output_path=output_path,
        analysis_logic_name="repeatmask_customlib",
        analysis_display_label="label",
        analysis_description="desc",
        program="RepeatMasker",
        program_version="4.1.5",
        program_parameters=None,
        source_provider="Ensembl",
        is_primary=True,
    )

    doc = json.loads(output_path.read_text(encoding="utf-8"))
    assert "program_parameters" not in doc["analysis"]
    assert doc["source"]["is_primary"] is True
