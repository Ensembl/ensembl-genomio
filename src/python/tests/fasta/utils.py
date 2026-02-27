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
"""Utility methods for testing of `ensembl.io.genomio.fasta modules`."""
from pathlib import Path

from Bio import SeqIO


def read_fasta(path: Path) -> dict[str, str]:
    with open(path, "r", encoding="utf-8") as fh:
        return {r.id: str(r.seq) for r in SeqIO.parse(fh, "fasta")}


def force_open_failure_for_suffix(suffix: str):
    """Monkeypatches `open` to raise an IOError when trying to open a file with the given suffix."""
    exception = OSError(f"Simulated open failure for files ending with '{suffix}'")
    real_open = open

    def patched_open(file, *args, **kwargs):
        if str(file).endswith(suffix):
            raise exception
        return real_open(file, *args, **kwargs)

    return patched_open
