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
"""Parse the Rfam hits TSV into GenomIO ncRNA records."""

__all__ = ["RfamConverter",
    "RfamParsedRow"]

import argparse
from dataclasses import dataclass
from pathlib import Path

from ensembl.io.genomio.features.convert_to_genomio_json.base import (
    ConverterOptions,
    FeatureConverter,
    ParseFeaturesResult,
    format_parse_errors,
    parse_token,
    register_converter,
    register_top_level_converter,
    validate_parsed_coordinates,
)
from ensembl.utils.archive import open_gz_file

RFAM_HITS_COLUMNS = 15

@register_top_level_converter
@register_converter
class RfamConverter(FeatureConverter):
    """Converter for the Rfam hits TSV."""

    analysis_logic_name = "cmscan_rfam_14.5"
    command = "rfam"
    ncrna_tool = "cmscan"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the Rfam subcommand parser."""
        rfam_parser = subparsers.add_parser(
            cls.command,
            help="Convert the Rfam hits TSV to GenomIO JSON.",
        )
        cls.add_common_arguments(rfam_parser)
        rfam_parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label="Rfam Models",
            analysis_description=(
                "Covariance models from <a href='https://rfam.xfam.org'>Rfam</a>, "
                "aligned to the genome with 'cmscan' from the "
                "<a href='http://eddylab.org/infernal'>Infernal</a> suite of programs."
            ),
            program="Infernal",
            source_provider="Rfam",
        )

    @classmethod
    def parse_features(cls, input_path: Path, _options: ConverterOptions | None = None) -> ParseFeaturesResult:
        """Parse the Rfam hits TSV."""
        return parse_output(input_path)


@dataclass(frozen=True)
class RfamParsedRow:
    """Parsed Rfam hits TSV row."""

    feature: dict[str, object]


def parse_row(input_path: Path, line: str) -> RfamParsedRow:
    """Parse one Rfam hits TSV row."""
    # The description is the final field and may contain spaces.
    columns = line.split(maxsplit=RFAM_HITS_COLUMNS - 1)
    if len(columns) != RFAM_HITS_COLUMNS:
        raise ValueError(
            f"Expected {RFAM_HITS_COLUMNS} columns in {input_path}, got {len(columns)}: line={line!r}"
        )

    seq_region, seq_start, seq_end, strand = columns[:4]
    truncation, gc, bias = columns[4:7]
    model_start, model_end = columns[7:9]
    score, evalue = columns[9:11]
    target_name, target_accession, biotype, description = columns[11:15]

    model_start_parse = parse_token(int, model_start, "model_start", line, input_path)
    model_end_parse = parse_token(int, model_end, "model_end", line, input_path)
    seq_start_parse = parse_token(int, seq_start, "seq_region_start", line, input_path)
    seq_end_parse = parse_token(int, seq_end, "seq_region_end", line, input_path)
    seq_region_start, seq_region_end = min(seq_start_parse, seq_end_parse), max(seq_start_parse, seq_end_parse)
    gc_parse = parse_token(float, gc, "gc", line, input_path)
    bias_parse = parse_token(float, bias, "bias", line, input_path)
    validate_parsed_coordinates(
        input_path,
        seq_region_start=seq_region_start,
        seq_region_end=seq_region_end,
        seq_region_strand=strand,
        repeat_start=model_start_parse,
        repeat_end=model_end_parse,
        line=line,
    )

    feature: dict[str, object] = {
        "seq_region": seq_region,
        "seq_region_start": seq_region_start,
        "seq_region_end": seq_region_end,
        "model_start": model_start_parse,
        "model_end": model_end_parse,
        "seq_region_strand": strand,
        "score": parse_token(float, score, "score", line, input_path),
        "evalue": parse_token(float, evalue, "evalue", line, input_path),
        "gc": gc_parse,
        "bias": bias_parse,
        "truncation": truncation,
        "biotype": biotype,
        "target_name": target_name,
        "description": description,
    }
    if target_accession != "-":
        feature["target_accession"] = target_accession
    return RfamParsedRow(feature)


def parse_output(input_path: Path) -> ParseFeaturesResult:
    """Parse all Rfam hits TSV rows, collating malformed-row errors."""
    features: list[dict[str, object]] = []
    errors: list[str] = []
    with open_gz_file(input_path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if line.lower().startswith("seqname"):
                continue
            try:
                parsed_row = parse_row(input_path, line)
            except ValueError as exc:
                errors.append(str(exc))
                continue

            features.append(parsed_row.feature)

    if errors:
        raise ValueError(format_parse_errors("Rfam hits TSV", input_path, errors))
    return features, {}
