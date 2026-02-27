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
"""Local directory-specific hook implementations for `ensembl.io.genomio.fasta` modules."""
import gzip
import pytest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def write_fasta(tmp_path):
    """Helper to create FASTA (optionally gzipped) from a list of (id, seq, desc)."""

    def _writer(relpath: str, record_info: list[tuple[str, str, str | None]], gz: bool = False) -> Path:
        p = tmp_path / relpath
        p.parent.mkdir(parents=True, exist_ok=True)

        records: list[SeqRecord] = []
        for record_id, seq, desc in record_info:
            records.append(SeqRecord(Seq(seq), id=record_id, description=desc or ""))

        opener = open
        mode = "w"
        if gz:
            opener = gzip.open
            mode = "wt"
            if p.suffix != ".gz":
                p = p.with_suffix(p.suffix + ".gz") if p.suffix else Path(str(p) + ".gz")
        with opener(p, mode, encoding="utf-8", newline="\n") as fh:
            SeqIO.write(records, fh, "fasta")

        return p

    return _writer
