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
"""Unit testing of `ensembl.io.genomio.utils.agp_utils` module."""

import pytest

from ensembl.io.genomio.utils import agp_utils


def test_parse_agp_ignores_comments_and_blank_lines(tmp_path):
    agp = tmp_path / "x.agp"
    agp.write_text(
        "# comment\n\n" "obj1\t1\t4\t1\tW\tpart1\t1\t4\t+\n",
        encoding="utf-8",
    )
    out = agp_utils.parse_agp(agp, allow_revcomp=False)
    assert "obj1" in out
    assert out["obj1"][0].part_id == "part1"


def test_parse_agp_rejects_short_line(tmp_path):
    agp = tmp_path / "bad.agp"
    agp.write_text("obj\t1\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected >= 9 columns"):
        agp_utils.parse_agp(agp, allow_revcomp=False)


def test_parse_agp_rejects_non_W_component(tmp_path):
    agp = tmp_path / "bad.agp"
    agp.write_text("obj\t1\t2\t1\tN\tgap\t1\t2\t+\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported AGP component type"):
        agp_utils.parse_agp(agp, allow_revcomp=False)


def test_parse_agp_rejects_minus_orientation_unless_allowed(tmp_path):
    agp = tmp_path / "bad.agp"
    agp.write_text("obj\t1\t2\t1\tW\tpart\t1\t2\t-\n", encoding="utf-8")
    with pytest.raises(ValueError, match="--allow-revcomp is not enabled"):
        agp_utils.parse_agp(agp, allow_revcomp=False)

    ok = agp_utils.parse_agp(agp, allow_revcomp=True)
    assert ok["obj"][0].orientation == "-"


def test_parse_agp_empty_file_raises(tmp_path):
    agp = tmp_path / "empty.agp"
    agp.write_text("# only comment\n", encoding="utf-8")
    with pytest.raises(ValueError, match="contained no component lines"):
        agp_utils.parse_agp(agp, allow_revcomp=False)
