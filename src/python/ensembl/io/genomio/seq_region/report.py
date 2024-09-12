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
"""Object for an INSDC assembly report to expose its data and metadata easily."""

__all__ = [
    "ReportRecord",
]

import csv
from os import PathLike
from pathlib import Path
import re
from typing import Tuple

from ensembl.utils.archive import open_gz_file


class ReportRecord:
    """Represent an assembly report file. Exposes 2 things:
    - Metadata as a dict from the comments.
    - A DictReader that yields all the seq_region lines of the report as dicts.
    """

    def __init__(self, report_path: Path) -> None:
        report_csv, metadata = self.report_to_csv(report_path)
        self.metadata = metadata
        self.reader = csv.DictReader(report_csv.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE)

    @staticmethod
    def report_to_csv(report_path: PathLike) -> Tuple[str, dict]:
        """Returns an assembly report as a CSV string.

        Args:
            report_path: Path to a seq_region file from INSDC/RefSeq.

        Returns:
            The data as a string in CSV format, and the head metadata as a dictionary.

        """
        with open_gz_file(report_path) as report:
            data = ""
            metadata = {}
            header_line = ""
            for line in report:
                if line.startswith("#"):
                    # Get metadata values if possible
                    match = re.search("# (.+?): (.+?)$", line)
                    if match:
                        metadata[match.group(1)] = match.group(2)
                    header_line = line
                else:
                    data += line

            if not header_line:
                raise ValueError("Missing header in report")
            data = header_line[2:].strip() + "\n" + data

            return data, metadata
