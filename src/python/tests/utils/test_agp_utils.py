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

from contextlib import nullcontext as does_not_raise
from pathlib import Path
import pytest
from pytest import param
from typing import ContextManager

from ensembl.io.genomio.utils import agp_utils


@pytest.mark.parametrize(
    "test_dir_name,agp_name,allow_revcomp,expectation",
    [
        param(
            "ignores_comments",
            "comment.agp",
            False,
            does_not_raise(("obj", "part", "+")),
            id="ignores_comments",
        ),
        param(
            "handles_orientation",
            "minus_strand.agp",
            True,
            does_not_raise(("obj", "part", "-")),
            id="handles_minus_strand_when_allowed",
        ),
        param(
            "handles_orientation",
            "minus_strand.agp",
            False,
            pytest.raises(ValueError, match=r"--allow-revcomp is not enabled"),
            id="rejects_minus_strand_when_not_allowed",
        ),
        param(
            "non_w_component",
            "repeat.agp",
            False,
            pytest.raises(ValueError, match=r"Unsupported AGP component type"),
            id="rejects_non_W_component",
        ),
        param(
            "empty_file",
            "empty.agp",
            False,
            pytest.raises(ValueError, match=r"contained no component lines"),
            id="rejects_empty_file",
        ),
        param(
            "truncated_line",
            "truncated.agp",
            False,
            pytest.raises(ValueError, match=r"expected >= 9 columns"),
            id="rejects_truncated_line",
        ),
    ],
)
def test_parse_agp(
    data_dir: Path,
    test_dir_name: str,
    agp_name: str,
    allow_revcomp: bool,
    expectation: ContextManager,
) -> None:
    """
    Tests the agp_utils.parse_agp() function.

    Args:
        data_dir: Module's test data directory fixture.
        test_dir_name: Name of the subdirectory within the test data directory that contains the
            AGP file being tested.
        agp_name: Name of the test AGP file.
        allow_revcomp: Boolean indicating whether minus strand entries are permitted.
        expectation: Context manager for the expected outcome of the test (exception or not).
    """
    test_dir = data_dir / test_dir_name
    agp_file = test_dir / agp_name

    with expectation as expected:
        out = agp_utils.parse_agp(agp_file, allow_revcomp)
        object_id, part_id, orientation = expected
        assert object_id in out
        assert out[object_id][0].part_id == part_id
        assert out[object_id][0].orientation == orientation


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
    "start,end,orientation,allow_revcomp,expectation",
    [
        param(
            1,
            1,
            "+",
            False,
            does_not_raise((100, 100)),
            id="single_bp_forward_strand",
        ),
        param(
            10,
            20,
            "+",
            False,
            does_not_raise((109, 119)),
            id="range_forward_strand",
        ),
        param(
            1,
            1,
            "-",
            True,
            does_not_raise((199, 199)),
            id="single_bp_reverse_strand",
        ),
        param(
            10,
            20,
            "-",
            True,
            does_not_raise((180, 190)),
            id="range_reverse_strand",
        ),
        param(
            5,
            4,
            "+",
            False,
            pytest.raises(ValueError, match=r"Range start > end"),
            id="rejects_start_after_end",
        ),
        param(
            0,
            1,
            "+",
            False,
            pytest.raises(ValueError, match=r"outside component span"),
            id="rejects_start_outside_component_span",
        ),
        param(
            1,
            101,
            "+",
            False,
            pytest.raises(ValueError, match=r"outside component span"),
            id="rejects_end_outside_component_span",
        ),
        param(
            1,
            1,
            "-",
            False,
            pytest.raises(ValueError, match=r"--allow-revcomp is not enabled"),
            id="rejects_minus_strand_when_not_allowed",
        ),
        param(
            1,
            10,
            "?",
            True,
            pytest.raises(ValueError, match=r"Invalid AGP orientation"),
            id="rejects_invalid_orientation",
        ),
    ],
)
def test_lift_range(
    start: int,
    end: int,
    orientation: str,
    allow_revcomp: bool,
    expectation: ContextManager,
) -> None:
    """
    Tests the agp_utils.lift_range() function.

    Args:
        start: Start position relative to component.
        end: End position relative to component.
        orientation: Orientation of the component in relation to the record (+/-).
        allow_revcomp: Boolean indicating whether minus strand entries are permitted.
        expectation: Context manager for the expected outcome of the test (exception or not).
    """
    with expectation as expected:
        part = agp_utils.AgpEntry(
            record="obj",
            record_start=100,
            record_end=199,
            part_number=1,
            part_id="comp",
            part_start=1,
            part_end=100,
            orientation=orientation,
        )

        obj, start_on_obj, end_on_obj = agp_utils.lift_range(part, start, end, allow_revcomp=allow_revcomp)
        assert obj == "obj"
        assert (start_on_obj, end_on_obj) == expected
