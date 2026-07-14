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
"""Unit testing of shared ``convert_to_genomio_json`` helpers."""

import argparse
from contextlib import nullcontext as does_not_raise
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from typing import Callable, ClassVar, ContextManager
from unittest.mock import Mock, patch

import pytest

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import base
from ensembl.io.genomio.features.convert_to_genomio_json.base import main, parse_args

from .helpers import sha256_key

def test_consensus_sha256_key_normalises_fields() -> None:
    """Test ``convert_to_genomio_json.Consensus.sha256_key()`` normalises whitespace and sequence case."""
    consensus = convert_to_genomio_json.Consensus(
        name=" AluY ",
        repeat_class=" SINE ",
        repeat_type=" Alu ",
        seq="ac gt\n",
    )

    assert consensus.sha256_key() == sha256_key("AluY", "SINE", "Alu", "ACGT")


@pytest.mark.parametrize(
    ("parser", "token", "field_name", "raw_line", "expectation"),
    [
        pytest.param(int, "42", "count", "42", does_not_raise(42), id="Valid integer"),
        pytest.param(
            int,
            "abc",
            "count",
            "abc",
            pytest.raises(ValueError, match=r"Invalid 'count' in input\.out: token='abc', line='abc'"),
            id="Invalid integer",
        ),
        pytest.param(float, "0.25", "score", "0.25", does_not_raise(0.25), id="Valid float"),
        pytest.param(
            float,
            "abc",
            "score",
            "abc",
            pytest.raises(ValueError, match=r"Invalid 'score' in input\.out: token='abc', line='abc'"),
            id="Invalid float",
        ),
    ],
)
def test_parse_token(
    parser: Callable[[str], int | float],
    token: str,
    field_name: str,
    raw_line: str,
    expectation: ContextManager,
) -> None:
    """Test ``base.parse_token()`` parses valid tokens and reports invalid token context.

    Args:
        parser: Parser callable to apply to the token.
        token: Raw token to parse.
        field_name: Field name used in error messages.
        raw_line: Original line used in error messages.
        expectation: Context manager for the expected result or exception.

    """
    with expectation as expected:
        assert base.parse_token(parser, token, field_name, raw_line, Path("input.out")) == expected


def test_format_parse_errors() -> None:
    """Test ``base.format_parse_errors()`` formats a counted bullet list."""
    output = base.format_parse_errors("dummy output", Path("input.dat"), ["first error", "second error"])
    assert output == "Found 2 errors while parsing dummy output in input.dat:\n- first error\n- second error"


def test_file_last_modified_time_returns_utc_isoformat(tmp_path: Path) -> None:
    """Test ``base.file_last_modified_time()`` returns a UTC ISO timestamp.

    Args:
        tmp_path: Temporary directory provided by pytest.

    """
    input_path = tmp_path / "input.out"
    input_path.write_text("content", encoding="utf-8")
    modified_time = datetime(2024, 1, 2, 3, 4, 5, tzinfo=timezone.utc).timestamp()

    os.utime(input_path, (modified_time, modified_time))

    assert base.file_last_modified_time(input_path) == "2024-01-02T03:04:05Z"


@pytest.mark.parametrize(
    ("seq_region_start", "seq_region_end", "repeat_start", "repeat_end", "expectation"),
    [
        pytest.param(1, 10, 2, 5, does_not_raise(), id="Valid coordinates"),
        pytest.param(
            0,
            10,
            1,
            5,
            pytest.raises(ValueError, match=r"Invalid seq_region coordinates"),
            id="Non-positive sequence region start",
        ),
        pytest.param(
            10,
            9,
            1,
            5,
            pytest.raises(ValueError, match=r"seq_region_end < seq_region_start"),
            id="Sequence region end before start",
        ),
        pytest.param(
            1,
            10,
            0,
            5,
            pytest.raises(ValueError, match=r"Invalid repeat coordinates"),
            id="Non-positive repeat start",
        ),
        pytest.param(
            1,
            10,
            5,
            4,
            pytest.raises(ValueError, match=r"repeat_end < repeat_start"),
            id="Repeat end before start",
        ),
    ],
)
def test_validate_parsed_coordinates(
    *,
    seq_region_start: int,
    seq_region_end: int,
    repeat_start: int,
    repeat_end: int,
    expectation: ContextManager,
) -> None:
    """Test ``base.validate_parsed_coordinates()`` correctly validates coordinates.

    Args:
        seq_region_start: Sequence region start coordinate.
        seq_region_end: Sequence region end coordinate.
        repeat_start: Repeat start coordinate.
        repeat_end: Repeat end coordinate.
        expectation: Context manager for the expected result or exception.

    """
    with expectation:
        base.validate_parsed_coordinates(
            Path("input.out"),
            seq_region_start=seq_region_start,
            seq_region_end=seq_region_end,
            repeat_start=repeat_start,
            repeat_end=repeat_end,
            line="raw line",
        )


class _DummyConverter(convert_to_genomio_json.FeatureConverter):
    """Converter used to test generic behavior without tool-specific fixtures."""

    analysis_logic_name = "dummy"
    analysis_display_label = "Dummy"
    analysis_description = "Dummy converter"
    command = "dummy"
    program = "dummy-program"
    consensuses_by_key: ClassVar[dict[str, convert_to_genomio_json.Consensus]] = {}

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add a generic dummy converter parser."""
        dummy_parser = subparsers.add_parser(cls.command)
        cls.add_common_arguments(dummy_parser)
        dummy_parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label=cls.analysis_display_label,
            analysis_description=cls.analysis_description,
            program=cls.program,
        )

    @classmethod
    def parse_features(
        cls,
        input_path: Path,
        _options: convert_to_genomio_json.ConverterOptions | None = None,
    ) -> convert_to_genomio_json.ParseFeaturesResult:
        """Return fixed parser output for document assembly tests."""
        return (
            [{"seq_region": input_path.name}],
            cls.consensuses_by_key,
        )


def test_converter_registration(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test generic registration helpers populate converter registries."""
    monkeypatch.setattr(base, "CONVERTERS_BY_LOGIC_NAME", {})
    monkeypatch.setattr(base, "TOP_LEVEL_CONVERTERS", [])

    assert base.register_converter(_DummyConverter) is _DummyConverter
    assert base.register_top_level_converter(_DummyConverter) is _DummyConverter
    assert base.register_converter(_DummyConverter) is _DummyConverter
    assert base.register_top_level_converter(_DummyConverter) is _DummyConverter
    assert {"dummy": _DummyConverter} == base.CONVERTERS_BY_LOGIC_NAME
    assert [_DummyConverter] == base.TOP_LEVEL_CONVERTERS


def test_duplicate_converter_logic_names_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test converter registration rejects conflicting classes for one logic name."""
    monkeypatch.setattr(base, "CONVERTERS_BY_LOGIC_NAME", {})

    class ConflictingConverter(_DummyConverter):
        """Converter with the same logic name as the dummy converter."""

    base.register_converter(_DummyConverter)

    with pytest.raises(ValueError, match=r"Converter logic name 'dummy' already registered"):
        base.register_converter(ConflictingConverter)


@pytest.mark.parametrize(
    ("analysis_logic_name", "converter_name"),
    [
        pytest.param(None, "NoLogicNameConverter", id="missing logic name"),
        pytest.param("", "BlankLogicNameConverter", id="blank logic name"),
    ],
)
def test_converters_without_logic_name_rejected(
    monkeypatch: pytest.MonkeyPatch,
    analysis_logic_name: str | None,
    converter_name: str,
) -> None:
    """Test converter registration requires a non-blank logic name."""
    monkeypatch.setattr(base, "CONVERTERS_BY_LOGIC_NAME", {})
    converter = type(converter_name, (_DummyConverter,), {"analysis_logic_name": analysis_logic_name})

    with pytest.raises(ValueError, match=rf"Converter {converter_name} has no analysis logic name"):
        base.register_converter(converter)


def test_duplicate_top_level_converter_commands_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test top-level converter registration rejects conflicting classes for one command."""
    monkeypatch.setattr(base, "TOP_LEVEL_CONVERTERS", [])

    class ConflictingTopLevelConverter(convert_to_genomio_json.FeatureConverter):
        """Top-level converter with the same command as the dummy converter."""

        command = "dummy"

        @classmethod
        def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
            """Add a generic conflicting parser."""

        @classmethod
        def parse_features(
            cls,
            input_path: Path,
            _options: convert_to_genomio_json.ConverterOptions | None = None,
        ) -> convert_to_genomio_json.ParseFeaturesResult:
            """Return fixed parser output for document assembly tests."""
            return (
                [{"seq_region": input_path.name}],
                {},
            )

    base.register_top_level_converter(_DummyConverter)

    with pytest.raises(ValueError, match=r"Top-level converter command 'dummy' already registered"):
        base.register_top_level_converter(ConflictingTopLevelConverter)


@pytest.mark.parametrize(
    ("command", "converter_name"),
    [
        pytest.param(None, "NoCommandTopLevelConverter", id="missing command"),
        pytest.param("", "BlankCommandTopLevelConverter", id="blank command"),
    ],
)
def test_top_level_converters_without_commands_rejected(
    monkeypatch: pytest.MonkeyPatch,
    command: str | None,
    converter_name: str,
) -> None:
    """Test top-level converter registration requires a non-blank CLI command."""
    monkeypatch.setattr(base, "TOP_LEVEL_CONVERTERS", [])
    converter = type(
        converter_name,
        (convert_to_genomio_json.FeatureConverter,),
        {
            "add_parser": classmethod(lambda _cls, _subparsers: None),
            "command": command,
        },
    )

    with pytest.raises(ValueError, match=rf"Top-level converter {converter_name} has no command"):
        base.register_top_level_converter(converter)


@pytest.mark.parametrize(
    "use_common_overrides",
    [
        pytest.param(False, id="Required common arguments only"),
        pytest.param(True, id="Optional common arguments overridden"),
    ],
)
def test_parse_args_common_arguments(
    convert_to_genomio_json_data_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    use_common_overrides: bool,
) -> None:
    """Test CLI arguments shared by converter subcommands."""
    program_version = "1.0"
    input_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    output_path = tmp_path / "out.json"
    monkeypatch.setattr(base, "TOP_LEVEL_CONVERTERS", [_DummyConverter])
    argv = [
        "dummy",
        "--input",
        str(input_path),
        "--output",
        str(output_path),
        "--program-version",
        program_version,
    ]

    if use_common_overrides:
        argv.extend(
            [
                "--program-parameters",
                "-lib foo",
                "--source-provider",
                "Custom",
                "--is-primary",
            ]
        )

    args = parse_args(argv)

    assert args.__class__.__name__ == "Namespace"
    assert args.input == input_path
    assert args.output == output_path
    assert args.program_version == program_version

    if use_common_overrides:
        assert args.source_provider == "Custom"
        assert args.is_primary is True
        assert args.program_parameters == "-lib foo"
    else:
        assert args.source_provider == "Ensembl"
        assert args.is_primary is False
        assert not hasattr(args, "program_parameters")


@patch("ensembl.io.genomio.features.convert_to_genomio_json.base.create_genomio_json")
def test_main_passes_common_config_fields(
    mock_create_genomio_json: Mock,
    convert_to_genomio_json_data_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test ``convert_to_genomio_json.main()`` passes common parsed config fields."""
    input_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    output_path = tmp_path / "out.json"
    monkeypatch.setattr(base, "TOP_LEVEL_CONVERTERS", [_DummyConverter])
    monkeypatch.setitem(base.CONVERTERS_BY_LOGIC_NAME, "dummy", _DummyConverter)

    main(
        [
            "dummy",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--program-version",
            "1.0",
            "--program-parameters",
            "params",
            "--source-provider",
            "Custom",
            "--is-primary",
        ]
    )

    mock_create_genomio_json.assert_called_once_with(
        config=convert_to_genomio_json.GenomioJsonConfig(
            input_path=input_path,
            output_path=output_path,
            analysis_logic_name="dummy",
            analysis_display_label="Dummy",
            analysis_description="Dummy converter",
            program="dummy-program",
            program_version="1.0",
            source_provider="Custom",
            is_primary=True,
            program_parameters="params",
        )
    )


@patch("ensembl.io.genomio.features.convert_to_genomio_json.base.create_genomio_json")
def test_main_reraises_exceptions(
    mock_create_genomio_json: Mock,
    convert_to_genomio_json_data_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test the ``convert_to_genomio_json.main()`` function reraises exceptions."""
    input_path = convert_to_genomio_json_data_dir / "create_json" / "basic.out"
    output_path = tmp_path / "out.json"
    monkeypatch.setattr(base, "TOP_LEVEL_CONVERTERS", [_DummyConverter])
    monkeypatch.setitem(base.CONVERTERS_BY_LOGIC_NAME, "dummy", _DummyConverter)

    mock_create_genomio_json.side_effect = RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        main(
            [
                "dummy",
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--program-version",
                "1.0",
            ]
        )


@pytest.mark.parametrize(
    ("consensuses_by_key", "expected_repeat_consensus"),
    [
        (
            {},
            None,
        ),
        (
            {
                "consensus-key": convert_to_genomio_json.Consensus(
                    name="Alu",
                    repeat_class="SINE",
                    repeat_type="Alu",
                    seq="ACGT",
                )
            },
            [
                {
                    "repeat_consensus_key": "consensus-key",
                    "repeat_name": "Alu",
                    "repeat_class": "SINE",
                    "repeat_type": "Alu",
                    "repeat_consensus": "ACGT",
                }
            ],
        ),
    ],
)
def test_create_genomio_json_assembles_generic_document(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    consensuses_by_key: dict[str, convert_to_genomio_json.Consensus],
    expected_repeat_consensus: list[dict[str, str]] | None,
) -> None:
    """Test JSON creation assembles parser output and generic metadata."""
    input_path = tmp_path / "input.out"
    input_path.write_text("parser input", encoding="utf-8")
    output_path = tmp_path / "out.json"
    monkeypatch.setitem(base.CONVERTERS_BY_LOGIC_NAME, "dummy", _DummyConverter)

    _DummyConverter.consensuses_by_key = consensuses_by_key
    convert_to_genomio_json.create_genomio_json(
        convert_to_genomio_json.GenomioJsonConfig(
            input_path=input_path,
            output_path=output_path,
            analysis_logic_name="dummy",
            analysis_display_label="Dummy",
            analysis_description="desc",
            program="dummy",
            program_version="1.0",
            source_provider="Ensembl",
            is_primary=False,
            program_parameters="-x",
        )
    )

    doc = json.loads(output_path.read_text(encoding="utf-8"))
    assert doc["analysis"]["logic_name"] == "dummy"
    assert doc["analysis"]["program_parameters"] == "-x"
    assert doc["source"] == {"source_provider": "Ensembl", "is_primary": False}
    assert doc["repeat_features"] == [{"seq_region": "input.out"}]
    if expected_repeat_consensus is None:
        assert "repeat_consensus" not in doc
    else:
        assert doc["repeat_consensus"] == expected_repeat_consensus


def test_create_genomio_json_rejects_unsupported_logic_name(
    convert_to_genomio_json_data_dir: Path, tmp_path: Path
) -> None:
    """Test JSON creation rejects unsupported analysis logic names."""
    with pytest.raises(ValueError, match=r"Unsupported analysis logic name"):
        convert_to_genomio_json.create_genomio_json(
            convert_to_genomio_json.GenomioJsonConfig(
                input_path=convert_to_genomio_json_data_dir / "create_json" / "basic.out",
                output_path=tmp_path / "out.json",
                analysis_logic_name="unsupported",
                analysis_display_label="label",
                analysis_description="desc",
                program="dummy",
                program_version="1.0",
                source_provider="Ensembl",
                is_primary=False,
            )
        )


def test_create_genomio_json_omits_program_parameters_when_none(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test that program parameters are omitted when no value is provided."""
    input_path = tmp_path / "input.out"
    input_path.write_text("parser input", encoding="utf-8")
    output_path = tmp_path / "out.json"
    monkeypatch.setitem(base.CONVERTERS_BY_LOGIC_NAME, "dummy", _DummyConverter)

    convert_to_genomio_json.create_genomio_json(
        convert_to_genomio_json.GenomioJsonConfig(
            input_path=input_path,
            output_path=output_path,
            analysis_logic_name="dummy",
            analysis_display_label="Dummy",
            analysis_description="desc",
            program="dummy",
            program_version="1.0",
            source_provider="Ensembl",
            is_primary=True,
        )
    )

    doc = json.loads(output_path.read_text(encoding="utf-8"))
    assert "program_parameters" not in doc["analysis"]
    assert doc["source"]["is_primary"] is True
