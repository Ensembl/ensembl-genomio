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

import pytest

from ensembl.io.genomio.features.convert_to_genomio_json import rfam

VALID_ROW = "NC_003076.8 4479446 4479542 - no 0.39 0.0 1 97 113.6 9.7e-27 MIR848 RF03308 pre_miRN MIR848 description"


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
    expectation: object,
    expected_feature: dict[str, object] | None,
) -> None:
    """Test that ``parse_row`` parses valid rows and rejects malformed rows."""
    with expectation:
        parsed_row = rfam.parse_row(Path("input.tsv"), line)

    if expected_feature is not None:
        assert parsed_row.feature == expected_feature


@pytest.mark.parametrize(
    ("expected_features", "expected_consensuses"),
    [
        pytest.param(
            [
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
            ],
            {},
            id="valid feature",
        )
    ],
)
def test_parse_output_success(
    tmp_path: Path,
    expected_features: list[dict[str, object]],
    expected_consensuses: dict[str, object],
) -> None:
    """Test that ``parse_output`` parses the header and one valid Rfam row."""
    input_file = tmp_path / "rfam_hits.tsv"
    input_file.write_text(
        "seqname\tstart\tend\tstrand\ttrunc\tgc\tbias\tmdl_from\tmdl_to\tscore\tevalue\tmodel_name\taccession\tbiotype\ttarget_description\n"
        f"{VALID_ROW}\n",
        encoding="utf-8",
    )
    features, consensuses_by_key = rfam.parse_output(input_file)

    assert features == expected_features
    assert consensuses_by_key == expected_consensuses

@pytest.mark.parametrize(
    ("invalid_rows", "error_pattern"),
    [
        pytest.param(
            [VALID_ROW.replace(" - ", " ? ", 1)],
            r"Found 1 errors while parsing Rfam hits TSV in .*:\n- Unexpected strand token.*",
            id="one parsing error",
        ),
        pytest.param(
            [
                VALID_ROW.replace(" - ", " ? ", 1),
                VALID_ROW.replace("4479446", "start", 1),
            ],
            (
                r"Found 2 errors while parsing Rfam hits TSV in .*:\n"
                r"- Unexpected strand token.*\n"
                r"- Invalid 'seq_region_start'.*"
            ),
            id="multiple parsing errors",
        ),
    ],
)
def test_parse_output_errors(tmp_path: Path, invalid_rows: list[str], error_pattern: str) -> None:
    """Test that ``parse_output`` aggregates row errors with ``format_parse_errors``."""
    input_file = tmp_path / "rfam_errors.tsv"
    input_file.write_text("\n".join(invalid_rows) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match=error_pattern):
        rfam.parse_output(input_file)
