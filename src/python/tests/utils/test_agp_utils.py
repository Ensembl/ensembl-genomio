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
import gzip

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


def test_parse_agp_decodes_bytes_from_gzipped_input(tmp_path):
    agp_gz = tmp_path / "x.agp.gz"

    payload = "# comment\nobj1\t1\t4\t1\tW\tpart1\t1\t4\t+\n"

    with gzip.open(agp_gz, "wt", encoding="utf-8") as fh:
        fh.write(payload)

    out = agp_utils.parse_agp(agp_gz, allow_revcomp=False)

    assert "obj1" in out
    assert len(out["obj1"]) == 1
    e = out["obj1"][0]
    assert e.part_id == "part1"
    assert e.orientation == "+"


def test_build_component_index_groups_by_part_id_and_sorts_within_component():
    e1 = agp_utils.AgpEntry(
        record="objB",
        record_start=50,
        record_end=60,
        part_number=2,
        part_id="comp1",
        part_start=1,
        part_end=11,
        orientation="+",
    )
    e2 = agp_utils.AgpEntry(
        record="objA",
        record_start=10,
        record_end=20,
        part_number=1,
        part_id="comp1",
        part_start=1,
        part_end=11,
        orientation="+",
    )
    e3 = agp_utils.AgpEntry(
        record="objA",
        record_start=21,
        record_end=30,
        part_number=2,
        part_id="comp1",
        part_start=1,
        part_end=10,
        orientation="+",
    )
    e4 = agp_utils.AgpEntry(
        record="objA",
        record_start=5,
        record_end=9,
        part_number=1,
        part_id="comp2",
        part_start=1,
        part_end=5,
        orientation="+",
    )

    agp_entries = {
        "objA": [e3, e2, e4],  # deliberately unsorted
        "objB": [e1],
    }

    idx = agp_utils.build_component_index(agp_entries)

    assert set(idx) == {"comp1", "comp2"}

    # comp2 only appears once
    assert idx["comp2"] == [e4]

    # comp1 should be sorted by (record, record_start, part_number):
    # objA start=10 pn=1, objA start=21 pn=2, objB start=50 pn=2
    assert idx["comp1"] == [e2, e3, e1]


def test_build_component_index_empty_input_returns_empty_dict():
    assert agp_utils.build_component_index({}) == {}


@pytest.mark.parametrize(
    "start,end,expected",
    [
        (1, 1, (100, 100)),
        (10, 20, (109, 119)),
        (100, 100, (199, 199)),
    ],
)
def test_lift_range_plus_orientation(start, end, expected):
    part = agp_utils.AgpEntry(
        record="obj",
        record_start=100,
        record_end=199,
        part_number=1,
        part_id="comp",
        part_start=1,
        part_end=100,
        orientation="+",
    )

    obj, s, e = agp_utils.lift_range(part, start, end, allow_revcomp=False)
    assert obj == "obj"
    assert (s, e) == expected


def test_lift_range_rejects_start_greater_than_end():
    part = agp_utils.AgpEntry(
        record="obj",
        record_start=1,
        record_end=10,
        part_number=1,
        part_id="comp",
        part_start=1,
        part_end=10,
        orientation="+",
    )
    with pytest.raises(ValueError, match=r"Range start > end"):
        agp_utils.lift_range(part, start=5, end=4, allow_revcomp=False)


def test_lift_range_rejects_range_outside_component_span():
    part = agp_utils.AgpEntry(
        record="obj",
        record_start=100,
        record_end=199,
        part_number=1,
        part_id="comp",
        part_start=10,
        part_end=109,
        orientation="+",
    )

    with pytest.raises(ValueError, match=r"outside component span"):
        agp_utils.lift_range(part, start=9, end=10, allow_revcomp=False)

    with pytest.raises(ValueError, match=r"outside component span"):
        agp_utils.lift_range(part, start=10, end=110, allow_revcomp=False)


def test_lift_range_minus_orientation_requires_allow_revcomp():
    part = agp_utils.AgpEntry(
        record="obj",
        record_start=100,
        record_end=199,
        part_number=1,
        part_id="comp",
        part_start=1,
        part_end=100,
        orientation="-",
    )

    with pytest.raises(ValueError, match=r"--allow-revcomp is not enabled"):
        agp_utils.lift_range(part, start=10, end=20, allow_revcomp=False)


@pytest.mark.parametrize(
    "start,end,expected",
    [
        (10, 20, (180, 190)),
        (1, 1, (199, 199)),
        (100, 100, (100, 100)),
    ],
)
def test_lift_range_minus_orientation_allowed(start, end, expected):
    part = agp_utils.AgpEntry(
        record="obj",
        record_start=100,
        record_end=199,
        part_number=1,
        part_id="comp",
        part_start=1,
        part_end=100,
        orientation="-",
    )

    obj, s, e = agp_utils.lift_range(part, start, end, allow_revcomp=True)
    assert obj == "obj"
    assert (s, e) == expected


def test_lift_range_rejects_invalid_orientation():
    part = agp_utils.AgpEntry(
        record="obj",
        record_start=1,
        record_end=10,
        part_number=1,
        part_id="comp",
        part_start=1,
        part_end=10,
        orientation="?",
    )
    with pytest.raises(ValueError, match=r"Invalid AGP orientation"):
        agp_utils.lift_range(part, start=1, end=2, allow_revcomp=True)


@pytest.mark.parametrize(
    "start,end",
    [
        (5, 4),  # start > end
        (0, 1),  # below part_start
        (1, 101),  # above part_end
    ],
)
def test_lift_range_invalid_ranges(start, end):
    part = agp_utils.AgpEntry(
        record="obj",
        record_start=1,
        record_end=100,
        part_number=1,
        part_id="comp",
        part_start=1,
        part_end=100,
        orientation="+",
    )

    with pytest.raises(ValueError):
        agp_utils.lift_range(part, start, end, allow_revcomp=False)
