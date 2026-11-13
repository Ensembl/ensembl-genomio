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
"""Unit testing of Rfam GenomIO JSON conversion helpers."""

from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import rfam

VALID_ROW = (
    "NC_003076.8 4479446 4479542 - no 0.39 0.0 1 97 113.6 9.7e-27 MIR848 RF03308 pre_miRN MIR848 description"
)


@pytest.mark.parametrize(
    ("line", "expectation", "expected_feature"),
    [
        pytest.param(
            VALID_ROW,
            does_not_raise(),
            {
                "seq_region": "NC_003076.8",
                "seq_region_start": 4479446,
                "seq_region_end": 4479542,
                "model_start": 1,
                "model_end": 97,
                "seq_region_strand": "-",
                "score": 113.6,
                "evalue": 9.7e-27,
                "gc": 0.39,
                "bias": 0.0,
                "truncation": "no",
                "biotype": "pre_miRN",
                "target_name": "MIR848",
                "target_accession": "RF03308",
                "description": "MIR848 description",
            },
            id="valid row",
        ),
        pytest.param(
            " ".join(VALID_ROW.split()[:14]),
            pytest.raises(ValueError, match="Expected 15 columns"),
            None,
            id="too many columns",
        ),
        pytest.param(
            "NC_003076.8 4479446 4479542",
            pytest.raises(ValueError, match="Expected 15 columns"),
            None,
            id="too few columns",
        ),
        pytest.param(
            VALID_ROW.replace(" - ", " ? ", 1),
            pytest.raises(ValueError, match="Unexpected strand token"),
            None,
            id="invalid strand",
        ),
        pytest.param(
            VALID_ROW.replace("4479446", "start", 1),
            pytest.raises(ValueError, match="Invalid 'seq_region_start'"),
            None,
            id="invalid sequence start",
        ),
        pytest.param(
            VALID_ROW.replace("4479446", "0", 1),
            pytest.raises(ValueError, match="Invalid seq_region coordinates"),
            None,
            id="non-positive sequence coordinate",
        ),
        pytest.param(
            VALID_ROW.replace(" 1 97 ", " first 97 ", 1),
            pytest.raises(ValueError, match="Invalid 'model_start'"),
            None,
            id="invalid model start",
        ),
        pytest.param(
            VALID_ROW.replace("113.6", "score", 1),
            pytest.raises(ValueError, match="Invalid 'score'"),
            None,
            id="invalid score",
        ),
        pytest.param(
            VALID_ROW.replace("9.7e-27", "evalue", 1),
            pytest.raises(ValueError, match="Invalid 'evalue'"),
            None,
            id="invalid evalue",
        ),
    ],
)
def test_parse_row(
    line: str,
    expectation: ContextManager[object],
    expected_feature: dict[str, object] | None,
) -> None:
    """Test that ``parse_row`` parses valid rows and rejects malformed rows.

    Args:
        line: Rfam row to parse.
        expectation: Expected success or exception context.
        expected_feature: Expected parsed feature, when parsing succeeds.

    Raises:
        ValueError: If an input row is malformed.

    """
    with expectation:
        parsed_row = rfam.parse_row(Path("input.tsv"), line)

    if expected_feature is not None:
        assert parsed_row.feature == expected_feature


def test_parse_output_success(
    tmp_path: Path,
    convert_to_genomio_json_data_dir: Path,
) -> None:
    """Test that ``parse_output`` parses the header and one valid Rfam row.

    Args:
        tmp_path: Pytest temporary directory fixture.
        convert_to_genomio_json_data_dir: Converter input fixture directory.

    """
    input_file = convert_to_genomio_json_data_dir / "rfam" / "success.tsv"
    config = convert_to_genomio_json.GenomioJsonConfig(
        input_path=input_file,
        output_path=tmp_path / "out.json",
        invalid_tool_output_path=tmp_path / "invalid.out",
        analysis_logic_name="rfam",
        analysis_display_label="Rfam",
        analysis_description="Rfam repeat finder",
        program="cmsearch",
        program_version="1.0",
        source_provider="Ensembl",
        is_primary=False,
        is_strict_mode=True,
        program_parameters=None,
        converter_options=convert_to_genomio_json.ConverterOptions(),
    )
    features, consensuses_by_key = rfam.parse_output(config)

    assert features == [
        {
            "seq_region": "NC_003076.8",
            "seq_region_start": 4479446,
            "seq_region_end": 4479542,
            "model_start": 1,
            "model_end": 97,
            "seq_region_strand": "-",
            "score": 113.6,
            "evalue": 9.7e-27,
            "gc": 0.39,
            "bias": 0.0,
            "truncation": "no",
            "biotype": "pre_miRN",
            "target_name": "MIR848",
            "target_accession": "RF03308",
            "description": "MIR848 description",
        }
    ]
    assert not consensuses_by_key


def test_rfam_converter_parse_features_delegates_to_parser(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test the Rfam converter delegates parsing using the complete config.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
        tmp_path: Pytest temporary directory fixture.

    """
    config = convert_to_genomio_json.GenomioJsonConfig(
        input_path=tmp_path / "input.tsv",
        output_path=tmp_path / "out.json",
        invalid_tool_output_path=tmp_path / "invalid.out",
        analysis_logic_name="rfam",
        analysis_display_label="Rfam",
        analysis_description="desc",
        program="cmsearch",
        program_version="1.0",
        source_provider="Rfam",
        is_primary=False,
        is_strict_mode=True,
    )
    expected: convert_to_genomio_json.ParseFeaturesResult = ([{"seq_region": "chr1"}], {})
    monkeypatch.setattr(rfam, "parse_output", lambda actual_config: expected)

    assert rfam.RfamConverter.parse_features(config) == expected


def test_parse_row_includes_target_accession_when_present() -> None:
    """Test that a non-placeholder target accession is retained."""
    row = VALID_ROW.replace("RF03308", "RF12345", 1)

    parsed = rfam.parse_row(Path("input.tsv"), row)

    assert parsed.feature["target_accession"] == "RF12345"


def test_parse_row_omits_placeholder_target_accession() -> None:
    """Test that the placeholder target accession is omitted."""
    row = VALID_ROW.replace("RF03308", "-", 1)

    parsed = rfam.parse_row(Path("input.tsv"), row)

    assert "target_accession" not in parsed.feature


@pytest.mark.parametrize(
    ("filename", "error_pattern"),
    [
        pytest.param(
            "errors_one.tsv",
            r"Found 1 errors while parsing Rfam hits TSV in .*:\n- Unexpected strand token.*",
            id="one parsing error",
        ),
        pytest.param(
            "errors_multiple.tsv",
            (
                r"Found 2 errors while parsing Rfam hits TSV in .*:\n"
                r"- Unexpected strand token.*\n"
                r"- Invalid 'seq_region_start'.*"
            ),
            id="multiple parsing errors",
        ),
    ],
)
def test_parse_output_errors(
    tmp_path: Path,
    convert_to_genomio_json_data_dir: Path,
    filename: str,
    error_pattern: str,
) -> None:
    """Test that ``parse_output`` aggregates row errors with ``format_parse_errors``.

    Args:
        tmp_path: Pytest temporary directory fixture.
        convert_to_genomio_json_data_dir: Converter input fixture directory.
        filename: Rfam fixture file containing malformed rows.
        error_pattern: Expected aggregated error message pattern.

    Raises:
        ValueError: If one or more rows fail parsing in strict mode.

    """
    input_file = convert_to_genomio_json_data_dir / "rfam" / filename

    config = convert_to_genomio_json.GenomioJsonConfig(
        input_path=input_file,
        output_path=tmp_path / "out.json",
        invalid_tool_output_path=tmp_path / "invalid.out",
        analysis_logic_name="red",
        analysis_display_label="RED",
        analysis_description="Red repeat finder",
        program="red",
        program_version="1.0",
        source_provider="Ensembl",
        is_primary=False,
        is_strict_mode=True,
        program_parameters=None,
        converter_options=convert_to_genomio_json.ConverterOptions(),
    )
    with pytest.raises(ValueError, match=error_pattern):
        rfam.parse_output(config)
