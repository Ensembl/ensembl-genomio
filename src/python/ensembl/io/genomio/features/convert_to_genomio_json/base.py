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
"""Generic framework for converting feature tool output to GenomIO JSON."""

__all__ = [
    "Consensus",
    "ConverterOptions",
    "FeatureConverter",
    "GenomioJsonConfig",
    "ParseFeaturesResult",
    "converters_by_logic_name",
    "create_genomio_json",
    "file_last_modified_time",
    "format_parse_errors",
    "parse_token",
    "register_converter",
    "register_top_level_converter",
    "top_level_converters",
    "validate_parsed_coordinates",
]

from abc import ABC, abstractmethod
import argparse
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
from inspect import isabstract
import json
import logging
from pathlib import Path
from typing import Any, Callable, TypeVar

import ensembl.io.genomio

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

T = TypeVar("T")
TopLevelConverterT = TypeVar("TopLevelConverterT", bound=type[Any])

# Shared data types


@dataclass(frozen=True)
class Consensus:
    """Repeat consensus record used for feature-to-consensus linking."""

    name: str
    repeat_class: str
    repeat_type: str
    seq: str

    def sha256_key(self) -> str:
        """Return a normalized SHA256 digest for this consensus record.

        The digest is computed from the consensus name, repeat class, repeat type,
        and normalized sequence content.

        Returns:
            SHA256 hex digest for the consensus record.

        """
        norm_name = self.name.strip()
        norm_class = self.repeat_class.strip()
        norm_type = self.repeat_type.strip()
        norm_seq = "".join(self.seq.split()).upper()
        payload = f"{norm_name}\t{norm_class}\t{norm_type}\t{norm_seq}".encode()
        return hashlib.sha256(payload).hexdigest()


ParseFeaturesResult = tuple[list[dict[str, object]], dict[str, Consensus]]


@dataclass(frozen=True)
class ConverterOptions:
    """Tool-specific options supplied to feature converters."""


@dataclass(frozen=True)
class GenomioJsonConfig:
    """Configuration for creating a GenomIO feature JSON document."""

    input_path: Path
    output_path: Path
    analysis_logic_name: str
    analysis_display_label: str
    analysis_description: str
    program: str
    program_version: str
    source_provider: str
    is_primary: bool
    program_parameters: str | None = None
    converter_options: ConverterOptions = field(default_factory=ConverterOptions)


# Converter interface and shared CLI arguments


class FeatureConverter(ABC):
    """Abstract class contract for tool-specific feature converters.

    Concrete converters declare their analysis metadata, add a command parser,
    and parse one tool output format. They are registered and invoked as
    classes rather than instantiated. See ``docs/user_guide/converter_modules.md``
    for instructions and an implementation template.
    """

    analysis_logic_name: str | None = None
    analysis_display_label: str | None = None
    analysis_description: str | None = None
    command: str | None = None
    program: str | None = None

    @classmethod
    @abstractmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add this converter's CLI parser."""
        raise NotImplementedError

    @staticmethod
    def add_common_arguments(subparser: ArgumentParser) -> None:
        """Add arguments shared by all supported analysis subcommands."""
        subparser.add_argument_src_path(
            "--input",
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
            help="Source provider for the features.",
        )
        subparser.add_argument(
            "--is-primary",
            action="store_true",
            help="Whether the source provider is the primary source for these features.",
        )
        subparser.add_log_arguments()

    @classmethod
    def options_from_args(cls, _args: argparse.Namespace) -> ConverterOptions:
        """Build converter-specific options from parsed command-line arguments."""
        return ConverterOptions()

    @classmethod
    @abstractmethod
    def parse_features(
        cls,
        input_path: Path,
        options: ConverterOptions | None = None,
    ) -> ParseFeaturesResult:
        """Parse features and consensus records from a tool output file."""
        raise NotImplementedError


def parse_token(parser: Callable[[str], T], token: str, field_name: str, raw_line: str, path: Path) -> T:
    """Parse a field from tool output into the specified class type.

    Args:
        parser: A callable that takes a string and returns a value of type T.
        token: Raw token to parse.
        field_name: Field name used in error messages.
        raw_line: Original input line.
        path: Input file path.

    Returns:
        The parsed value.

    Raises:
        ValueError: If the token cannot be cast to the specified type.

    """
    try:
        return parser(token)
    except ValueError as exc:
        raise ValueError(f"Invalid {field_name!r} in {path}: token={token!r}, line={raw_line!r}") from exc


def format_parse_errors(parser_name: str, input_path: Path, errors: list[str]) -> str:
    """Format multiple parser errors into a single exception message."""
    return f"Found {len(errors)} errors while parsing {parser_name} in {input_path}:\n" + "\n".join(
        f"- {error}" for error in errors
    )


def file_last_modified_time(file_path: Path) -> str:
    """Return the last modified time of the given file."""
    return (
        datetime.fromtimestamp(
            file_path.stat().st_mtime,
            tz=timezone.utc,
        )
        .isoformat()
        .replace("+00:00", "Z")
    )


def validate_parsed_coordinates(
    input_path: Path,
    *,
    seq_region_start: int,
    seq_region_end: int,
    repeat_start: int,
    repeat_end: int,
    line: str,
) -> None:
    """Validate parsed coordinate values for a feature.

    Args:
        input_path: Input file path used.
        seq_region_start: Start coordinate on the sequence region.
        seq_region_end: End coordinate on the sequence region.
        repeat_start: Start coordinate on the repeat consensus.
        repeat_end: End coordinate on the repeat consensus.
        line: Original input line for error reporting.

    Raises:
        ValueError: If sequence region or repeat coordinate values are invalid (i.e. negative, zero,
        or end < start).

    """
    if seq_region_start < 1 or seq_region_end < 1:
        raise ValueError(
            f"Invalid seq_region coordinates in {input_path}: "
            f"start={seq_region_start}, end={seq_region_end}, line={line!r}"
        )
    if seq_region_end < seq_region_start:
        raise ValueError(
            f"seq_region_end < seq_region_start in {input_path}: "
            f"start={seq_region_start}, end={seq_region_end}, line={line!r}"
        )

    if repeat_start < 1 or repeat_end < 1:
        raise ValueError(
            f"Invalid repeat coordinates in {input_path}: "
            f"repeat_start={repeat_start}, repeat_end={repeat_end}, line={line!r}"
        )
    if repeat_end < repeat_start:
        raise ValueError(
            f"repeat_end < repeat_start in {input_path}: "
            f"repeat_start={repeat_start}, repeat_end={repeat_end}, line={line!r}"
        )


# Converter registry

# Register every concrete converter class that parses a single analysis output,
# keyed by the analysis.logic_name written to the JSON document.
converters_by_logic_name: dict[str, type[FeatureConverter]] = {}

# Register only converters that should appear as first-level CLI tool commands
# under the main parser. A top-level converter may parse records itself, or it
# may only create nested mode-specific subcommands.
top_level_converters: list[type[Any]] = []


def _converter_signature(converter: type[Any]) -> tuple[object, ...]:
    """Return the registration-relevant attributes of a converter class.

    The signature defines when independently created converter classes may
    replace one another in a registry. It accepts composite top-level
    converters, which need only define ``command`` and ``add_parser``; absent
    FeatureConverter-specific attributes are represented by ``None``.

    Args:
        converter: Converter class to describe.

    Returns:
        Tuple containing converter metadata and implementation functions that
        affect command-line parsing or feature conversion.

    """

    def implementation(name: str) -> object:
        method = getattr(converter, name, None)
        return getattr(method, "__func__", method)

    return (
        getattr(converter, "analysis_logic_name", None),
        getattr(converter, "analysis_display_label", None),
        getattr(converter, "analysis_description", None),
        getattr(converter, "command", None),
        getattr(converter, "program", None),
        implementation("add_parser"),
        implementation("add_common_arguments"),
        implementation("options_from_args"),
        implementation("parse_features"),
    )


def register_converter(converter: type[FeatureConverter]) -> type[FeatureConverter]:
    """Register a converter by analysis logic name."""
    if isabstract(converter):
        raise ValueError(f"Cannot register abstract converter {converter.__name__}")

    if not converter.analysis_logic_name:
        raise ValueError(f"Converter {converter.__name__} has no analysis logic name")
    registered_converter = converters_by_logic_name.get(converter.analysis_logic_name)

    if registered_converter is not None and _converter_signature(
        registered_converter
    ) != _converter_signature(converter):
        raise ValueError(
            f"Converter logic name {converter.analysis_logic_name!r} already registered "
            f"by {registered_converter.__name__}"
        )
    converters_by_logic_name[converter.analysis_logic_name] = converter
    return converter


def register_top_level_converter(converter: TopLevelConverterT) -> TopLevelConverterT:
    """Register a converter that should appear as a top-level CLI command."""
    converter_command = getattr(converter, "command", None)
    if not converter_command:
        raise ValueError(f"Top-level converter {converter.__name__} has no command")
    if isabstract(converter):
        raise ValueError(f"Cannot register abstract top-level converter {converter.__name__}")
    for index, registered_converter in enumerate(top_level_converters):
        registered_command = getattr(registered_converter, "command", None)
        if registered_command == converter_command:
            if _converter_signature(registered_converter) != _converter_signature(converter):
                raise ValueError(
                    f"Top-level converter command {converter_command!r} already registered "
                    f"by {registered_converter.__name__}"
                )
            top_level_converters[index] = converter
            return converter
    top_level_converters.append(converter)
    return converter


# Argument parsing


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
    for converter in top_level_converters:
        converter.add_parser(subparsers)

    args = parser.parse_args(argv)
    init_logging_with_args(args)
    return args


# JSON document creation


def create_genomio_json(config: GenomioJsonConfig) -> None:
    """Create a GenomIO JSON document from feature identification tool output.

    Args:
        config: Configuration describing the input, output, analysis metadata, source metadata,
            and converter-specific options.

    Raises:
        ValueError: If an unsupported analysis logic name is provided.

    """
    try:
        converter = converters_by_logic_name[config.analysis_logic_name]
    except KeyError:
        raise ValueError(f"Unsupported analysis logic name: {config.analysis_logic_name}") from None
    features, consensuses_by_key = converter.parse_features(config.input_path, config.converter_options)

    analysis: dict[str, str] = {
        "run_date": file_last_modified_time(config.input_path),
        "logic_name": config.analysis_logic_name,
        "display_label": config.analysis_display_label,
        "description": config.analysis_description,
        "program": config.program,
        "program_version": config.program_version,
    }
    if config.program_parameters is not None:
        analysis["program_parameters"] = config.program_parameters

    json_doc: dict[str, object] = {
        "analysis": analysis,
        "source": {
            "source_provider": config.source_provider,
            "is_primary": config.is_primary,
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

    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    config.output_path.write_text(json.dumps(json_doc, indent=2) + "\n", encoding="utf-8")


# CLI entry point


def main(argv: list[str] | None = None) -> None:
    """Run the JSON conversion command-line entry point.

    Args:
        argv: Optional list of command-line arguments. If `None`, arguments are taken from ``sys.argv``.

    """
    args = parse_args(argv)
    try:
        converter = converters_by_logic_name[args.analysis_logic_name]
        create_genomio_json(
            config=GenomioJsonConfig(
                input_path=args.input,
                output_path=args.output,
                analysis_logic_name=args.analysis_logic_name,
                analysis_display_label=args.analysis_display_label,
                analysis_description=args.analysis_description,
                program=args.program,
                program_version=args.program_version,
                source_provider=args.source_provider,
                is_primary=args.is_primary,
                program_parameters=getattr(args, "program_parameters", None),
                converter_options=converter.options_from_args(args),
            )
        )
    except Exception:
        logging.exception(f"Error processing file {args.input}")
        raise
