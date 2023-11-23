#!/usr/bin/env python
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Checks if every flat file of a given format in the repository is valid."""

import argparse
import json
from pathlib import Path
import sys

import yaml


# Do not show traceback for an easier error readability
sys.tracebacklimit = 0


_ROOT_PATH = Path(__file__).absolute().parents[2]
_FORMAT_CONFIG = {
    "json": {
        "function": json.load,
        "kwargs": {},
        "exception": json.JSONDecodeError,
        "extensions": set(["json"]),
    },
    "yaml": {
        "function": yaml.load,
        "kwargs": {
            "Loader": yaml.Loader,
        },
        "exception": yaml.YAMLError,
        "extensions": set(["yml", "yaml"]),
    },
}


# Add GitLab's "!reference" and mkdocs's "!ENV" tags to yaml.Loader class
def yaml_tag_constructor(loader, node):
    """Default constructor to get the value of new YAML tags."""
    seq = loader.construct_sequence(node)
    return seq


yaml.Loader.add_constructor("!reference", yaml_tag_constructor)
yaml.Loader.add_constructor("!ENV", yaml_tag_constructor)


def check_flatfile(file_format: str) -> None:
    """Checks if every flat file in the repository for the given format is valid.

    Args:
        file_format: Flat file format to validate.

    Raises:
        RuntimeError: If at least one file is not valid.

    """
    config = _FORMAT_CONFIG[file_format]
    file_extensions = config["extensions"]
    report_files = []
    for file_path in _ROOT_PATH.rglob("*.*"):
        if file_path.is_file() and (file_path.suffix[1:] in file_extensions):
            with file_path.open() as fp:
                try:
                    config["function"](fp, **config["kwargs"])
                except config["exception"] as exc:
                    print(exc, file=sys.stderr)
                    report_files.append(str(file_path))
    if report_files:
        report = "\n".join(report_files)
        raise RuntimeError(f"{len(report_files)} flat files are not well formatted\n\n{report}")
    else:
        print(f"Every {file_format.upper()} file in the repository is well formatted")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validates all the repository flat files for the selected format."
    )
    parser.add_argument("file_format", choices=_FORMAT_CONFIG.keys(), help="Flat file format to validate")
    args = parser.parse_args()
    check_flatfile(args.file_format)
