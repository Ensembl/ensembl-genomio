# Adding a feature converter

Feature converters transform a tool-specific feature file into a GenomIO JSON document. The shared command-line interface, validation, registration, and document assembly code lives in `ensembl.io.genomio.features.convert_to_genomio_json.base`; each supported tool has its own module in `src/python/ensembl/io/genomio/features/convert_to_genomio_json/`.

## Implement a concrete converter

Create a module named for the tool, for example `my_tool.py`. Define one `FeatureConverter` subclass for each analysis that can produce a document. A concrete converter must provide:

- `analysis_logic_name`, `analysis_display_label`, `analysis_description`, `command`, and `program` class attributes;
- `add_parser()`, which creates the command parser and sets the analysis metadata defaults;
- `parse_features()`, which reads the tool output and returns a `ParseFeaturesResult`;
- `options_from_args()` only when the tool has options beyond the shared arguments.

Use `add_common_arguments()` on leaf analysis parsers. It adds the input, output, program metadata, source metadata, and logging arguments shared by all feature converters.

```python
import argparse
from pathlib import Path

from ensembl.io.genomio.features.convert_to_genomio_json.base import (
    ConverterOptions,
    FeatureConverter,
    ParseFeaturesResult,
    register_converter,
    register_top_level_converter,
)


@register_top_level_converter
@register_converter
class MyToolConverter(FeatureConverter):
    """Convert MyTool output to GenomIO JSON."""

    analysis_logic_name = "my_tool"
    analysis_display_label = "MyTool features"
    analysis_description = "Features reported by MyTool."
    command = "my-tool"
    program = "my-tool"

    @classmethod
    def add_parser(cls, subparsers: argparse._SubParsersAction) -> None:
        """Add the MyTool subcommand parser."""
        parser = subparsers.add_parser(cls.command, help="Convert MyTool output to GenomIO JSON.")
        cls.add_common_arguments(parser)
        parser.set_defaults(
            analysis_logic_name=cls.analysis_logic_name,
            analysis_display_label=cls.analysis_display_label,
            analysis_description=cls.analysis_description,
            program=cls.program,
        )

    @classmethod
    def parse_features(
        cls,
        input_path: Path,
        _options: ConverterOptions | None = None,
    ) -> ParseFeaturesResult:
        """Parse MyTool output."""
        return parse_my_tool_output(input_path)
```

`FeatureConverter` is abstract. Do not instantiate it or a concrete converter: converters are registered and invoked as classes.

## Registering converters

Use `@register_converter` for every concrete analysis converter. It registers the class by `analysis_logic_name`, which is used to select the parser when assembling JSON.

Also use `@register_top_level_converter` when the converter should appear as a first-level CLI command. A converter can be registered for parsing without being a top-level command.

Some tools have several modes beneath one top-level command. For these, create a composite top-level class that defines `command` and `add_parser()` but does not inherit `FeatureConverter`; `RepeatMaskerConverter` is the existing example. Register the concrete mode converters with `@register_converter`, and register only the composite class with `@register_top_level_converter`.

The registries accept a later class with the same registration signature as a replacement. A different signature using an existing `analysis_logic_name` or top-level `command` is rejected.

## Test a new converter

Add tool-specific tests alongside the converter:

- `src/python/tests/features/test_convert_to_genomio_json_<tool>.py` for parser arguments, parsed features, options, and malformed input;
- fixture data under `src/python/tests/features/test_convert_to_genomio_json/<tool>/` when representative tool output is needed.

Keep generic registry and shared-CLI behavior in `test_convert_to_genomio_json_base.py`. Run the focused test module while developing, then run the complete suite before submitting the change.

```bash
pytest src/python/tests/features/test_convert_to_genomio_json_my_tool.py
pytest
```
