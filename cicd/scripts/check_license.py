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
"""Checks if one of the requirements of the license is met by the repository."""

import argparse
import datetime
from itertools import zip_longest
from os import PathLike
from pathlib import Path
import re
import sys


# Do not show traceback for an easier error readability
sys.tracebacklimit = 0


_SCRIPTS_DIR = Path(__file__).absolute().parent
_TEMPLATES_DIR = _SCRIPTS_DIR.parent / "templates"
_ROOT_PATH = _SCRIPTS_DIR.parents[1]
# Set of file suffixes for which we expect to find the license header
_SUFFIXES_WITH_HEADER = set(
    ["pl", "pm", "py", "nf", "config", "mysql", "pgsql", "sql", "sqlite", "bash", "sh", "toml", "yml"]
)
_EXCLUDE_DIRS = ["data/test", "src/python/tests/data"]


def check_notice(notice_template: PathLike) -> None:
    """Checks if the NOTICE file is correct and the copyright year is up-to-date.

    Args:
        notice_template: Path to notice template file.

    Raises:
        RuntimeError: If the NOTICE file has incorrect format or the copyright year is not correct.

    """
    notice_file = _ROOT_PATH / "NOTICE"
    year = datetime.date.today().year
    report = ""
    with Path(notice_template).open() as tpl, notice_file.open() as fh:
        # Add dummy values if one of the files turns out to be shorter than the other
        for tpl_line, fh_line in zip_longest(tpl, fh, fillvalue=""):
            tpl_line = tpl_line.replace("<current_year>", f"{year}")
            tpl_line = tpl_line.rstrip("\n")
            fh_line = fh_line.rstrip("\n")
            if tpl_line != fh_line:
                report += f"> Expected: '{tpl_line}'\n  Found:    '{fh_line}'\n"
    if report:
        raise RuntimeError(f"Incorrect NOTICE file format or copyright year\n\n{report}")
    else:
        print("NOTICE file is correct and the copyright year is up-to-date")


def check_header(header_template: PathLike) -> None:
    """Checks if every code file in the repository has a valid license header.

    Args:
        header_template: Path to license header template file.

    Raises:
        RuntimeError: If at least one file has missing or incorrect license header.

    """
    template = Path(header_template).read_text()
    # Escape symbols that may be interpreted as part of a regex otherwise
    template = template.replace(".", "\.").replace("(", "\(").replace(")", "\)")
    # Allow comment symbols and additional spaces before each header line, and single newline instead of two
    template = "[#/-]*\s*" + template.replace("\n\n", "(\n)+").replace("\n", "\n[#/-]*\s*")
    # Compile template as regex to improve performance
    prog = re.compile(rf"{template}")
    report_files = []
    for file_path in _ROOT_PATH.rglob("*.*"):
        if file_path.is_file() and (file_path.suffix[1:] in _SUFFIXES_WITH_HEADER):
            if any([file_path.match(f"{_ROOT_PATH}/{x}/*") for x in _EXCLUDE_DIRS]):
                # Do not check any files that belong to one of the directories to exclude
                continue
            if not prog.search(file_path.read_text()):
                report_files.append(str(file_path))
    if report_files:
        report = "\n".join(report_files)
        raise RuntimeError(
            f"{len(report_files)} code files have missing or incorrect license header\n\n{report}"
        )
    else:
        print("Every code file in the repository has a valid license header")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Checks if the NOTICE file is correct or if every code file in the repository has a valid "
            "license header."
        )
    )
    # Create separate subparsers for each case "notice" and "header"
    subparsers = parser.add_subparsers(title="License aspect to check", required=True, dest="{notice,header}")
    parser_notice = subparsers.add_parser("notice", help="Check NOTICE file format and copyright year")
    parser_notice.add_argument(
        "--template",
        type=Path,
        required=False,
        default=_TEMPLATES_DIR / "notice.tpl",
        help="Path to notice template file",
    )
    parser_notice.set_defaults(check_function=check_notice)
    parser_header = subparsers.add_parser("header", help="Check if code files have proper license header")
    parser_header.add_argument(
        "--template",
        type=Path,
        required=False,
        default=_TEMPLATES_DIR / "license_header.tpl",
        help="Path to license header template file",
    )
    parser_header.set_defaults(check_function=check_header)
    args = parser.parse_args()
    # Run the corresponding function depending on the type of check selected by the user
    args.check_function(args.template)
