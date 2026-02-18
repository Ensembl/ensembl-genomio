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
"""Unit testing of `ensembl.io.genomio.repeats.combine_json` module."""

import gzip
import json
from pathlib import Path
import re

import pytest

from ensembl.io.genomio.repeats import combine_json  # <-- adjust if needed
from .._helpers import write_manifest


CHUNK_RE = re.compile(r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$")


@pytest.fixture()
def schema_validator_calls(monkeypatch: pytest.MonkeyPatch):
    """
    Records calls made to schema_validator and replaces it with a no-op.

    Returned value is a list of (args, kwargs) tuples for inspection in tests.
    """
    calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def _schema_validator_recorder(*args: object, **kwargs: object) -> None:
        calls.append((args, kwargs))

    monkeypatch.setattr(combine_json, "schema_validator", _schema_validator_recorder)
    return calls


def _write_json(path: Path, obj: object) -> Path:
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    return path


def _write_json_gz(path: Path, obj: object) -> Path:
    data = (json.dumps(obj, indent=2) + "\n").encode("utf-8")
    with gzip.open(path, "wb") as fh:
        fh.write(data)
    return path


def _analysis(name: str = "rm") -> dict[str, object]:
    return {"logic_name": name}


def _source(provider: str = "prov") -> dict[str, object]:
    return {"source_provider": provider, "is_primary": True}


def _consensus(rn: str = "Alu", rc: str = "SINE", rt: str = "Alu") -> combine_json.RepeatConsensus:
    return {"repeat_name": rn, "repeat_class": rc, "repeat_type": rt}


def _feature(
    seq_region: str,
    s: int,
    e: int,
    strand: int = 1,
    consensus: combine_json.RepeatConsensus | None = None,
) -> combine_json.RepeatFeature:
    feat: combine_json.RepeatFeature = {
        "seq_region": seq_region,
        "seq_region_start": s,
        "seq_region_end": e,
        "seq_region_strand": strand,
        "repeat_start": 1,
        "repeat_end": 10,
    }
    if consensus is not None:
        feat["repeat_consensus"] = consensus
    return feat


@pytest.mark.parametrize("gz", [False, True])
def test_load_json_document_accepts_object(tmp_path: Path, gz: bool):
    doc = {"analysis": _analysis(), "source": _source(), "repeat_features": []}
    p = tmp_path / ("in.json.gz" if gz else "in.json")
    (_write_json_gz if gz else _write_json)(p, doc)

    got = combine_json._load_json_document(p)
    assert isinstance(got, dict)
    assert "analysis" in got


def test_load_json_document_rejects_non_object(tmp_path: Path):
    p = tmp_path / "bad.json"
    _write_json(p, [1, 2, 3])

    with pytest.raises(ValueError, match=r"Top-level JSON .* must be an object"):
        combine_json._load_json_document(p)


@pytest.mark.parametrize("gz", [False, True])
def test_schema_validate_file_calls_validator(schema_validator_calls, tmp_path: Path, gz: bool):
    doc = {"analysis": _analysis(), "source": _source(), "repeat_features": []}
    p = tmp_path / ("v.json.gz" if gz else "v.json")
    (_write_json_gz if gz else _write_json)(p, doc)

    # Force the "tmp exists" branch for gz: v.json already exists.
    if gz:
        tmp_path.joinpath("v.json").write_text("preexisting", encoding="utf-8")

    combine_json._schema_validate_file(p)

    assert schema_validator_calls, "schema_validator was not called"
    args, kwargs = schema_validator_calls[-1]
    assert kwargs.get("json_schema") == "repeat"

    validated_path = kwargs.get("json_file") or (args[0] if args else None)
    assert isinstance(validated_path, Path)


def test_get_agp_entry_for_range_raises_when_no_span_matches(tmp_path: Path):
    AgpEntry = combine_json.AgpEntry
    parts = [
        AgpEntry("obj", 1, 100, 1, "comp", 1, 100, "+"),
    ]
    with pytest.raises(ValueError, match=r"does not fit within any AGP span"):
        combine_json._get_agp_entry_for_range(parts, 200, 210, component_id="comp", path=tmp_path / "x.agp")


def test_get_agp_entry_for_range_raises_when_ambiguous(tmp_path: Path):
    AgpEntry = combine_json.AgpEntry
    parts = [
        AgpEntry("obj", 1, 100, 1, "comp", 1, 100, "+"),
        AgpEntry("obj", 1, 100, 2, "comp", 50, 150, "+"),  # overlaps
    ]
    with pytest.raises(ValueError, match=r"Ambiguous AGP mapping"):
        combine_json._get_agp_entry_for_range(parts, 60, 70, component_id="comp", path=tmp_path / "x.agp")


def test_merge_repeat_consensus_rejects_conflicting_duplicates(tmp_path: Path):
    rc1: dict[str, combine_json.JsonValue] = {
        "repeat_name": "Alu",
        "repeat_class": "SINE",
        "repeat_type": "Alu",
    }
    rc2: dict[str, combine_json.JsonValue] = {
        "repeat_name": "Alu",
        "repeat_class": "SINE",
        "repeat_type": "Alu",
        "repeat_consensus": "ACGT",
    }

    combined: dict[tuple[str, str, str], combine_json.RepeatConsensus] = {}
    combine_json._merge_repeat_consensus(combined, [rc1], source_path=tmp_path / "a.json")

    with pytest.raises(ValueError, match=r"Conflicting repeat_consensus"):
        combine_json._merge_repeat_consensus(combined, [rc2], source_path=tmp_path / "b.json")


@pytest.mark.parametrize(
    "name,s,e,expected_region,expected_start,expected_end",
    [
        ("chr1_chunk_start_11", 1, 5, "chr1", 11, 15),
        ("chr1", 1, 5, "chr1", 1, 5),
    ],
)
def test_lift_repeat_feature_header_driven(
    tmp_path: Path,
    name: str,
    s: int,
    e: int,
    expected_region: str,
    expected_start: int,
    expected_end: int,
):
    feat = _feature(seq_region=name, s=s, e=e)

    out = combine_json._lift_repeat_feature(
        feat,
        chunk_re=CHUNK_RE,
        agp_by_component=None,
        allow_revcomp=False,
        source_path=tmp_path / "x.json",
    )

    assert out["seq_region"] == expected_region
    assert out["seq_region_start"] == expected_start
    assert out["seq_region_end"] == expected_end


def test_lift_repeat_feature_agp_driven_forward_and_reverse_strand(tmp_path: Path):
    AgpEntry = combine_json.AgpEntry

    forward_strand_entry = AgpEntry(
        record="objP",
        record_start=100,
        record_end=199,
        part_number=1,
        part_id="compP",
        part_start=1,
        part_end=100,
        orientation="+",
    )
    reverse_strand_entry = AgpEntry(
        record="objM",
        record_start=100,
        record_end=199,
        part_number=1,
        part_id="compM",
        part_start=1,
        part_end=100,
        orientation="-",
    )

    forward_strand_feature = _feature(seq_region="compP", s=10, e=20, strand=1)
    forward_strand_result = combine_json._lift_repeat_feature(
        forward_strand_feature,
        chunk_re=CHUNK_RE,
        agp_by_component={"compP": [forward_strand_entry]},
        allow_revcomp=False,
        source_path=tmp_path / "x.json",
    )
    assert forward_strand_result["seq_region"] == "objP"
    assert forward_strand_result["seq_region_start"] == 109
    assert forward_strand_result["seq_region_end"] == 119
    assert forward_strand_result["seq_region_strand"] == 1

    reverse_strand_feature = _feature(seq_region="compM", s=10, e=20, strand=1)
    reverse_strand_result = combine_json._lift_repeat_feature(
        reverse_strand_feature,
        chunk_re=CHUNK_RE,
        agp_by_component={"compM": [reverse_strand_entry]},
        allow_revcomp=True,
        source_path=tmp_path / "x.json",
    )
    assert reverse_strand_result["seq_region"] == "objM"
    assert reverse_strand_result["seq_region_strand"] == -1


def test_combine_repeat_json_header_driven_merge_and_coord_liftover(schema_validator_calls, tmp_path: Path):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _consensus("Alu", "SINE", "Alu")

    d1 = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_feature(seq_region="chr1_chunk_start_1", s=1, e=3, consensus=rc)],
    }
    d2 = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_feature(seq_region="chr1_chunk_start_4", s=1, e=2, consensus=rc)],
    }

    p1 = _write_json(tmp_path / "a.json", d1)
    p2 = _write_json(tmp_path / "b.json", d2)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p1), str(p2)])
    out_json = tmp_path / "out.json"

    combine_json.combine_repeat_json(
        json_manifest=manifest,
        out_json=out_json,
        chunk_re=CHUNK_RE,
        validate_inputs=True,
        validate_output=True,
    )

    out = json.loads(out_json.read_text(encoding="utf-8"))
    assert out["analysis"] == analysis
    assert out["source"] == source
    assert len(out["repeat_consensus"]) == 1
    assert len(out["repeat_features"]) == 2
    assert out["repeat_features"][0]["seq_region"] == "chr1"
    assert out["repeat_features"][1]["seq_region_start"] == 4

    assert schema_validator_calls


@pytest.mark.parametrize(
    "orientation,allow_revcomp,expected_strand,expected_start,expected_end",
    [
        ("+", False, 1, 109, 119),
        ("-", True, -1, 180, 190),
    ],
)
def test_combine_repeat_json_agp_driven_lifts_coordinates_and_strand(
    tmp_path: Path,
    orientation: str,
    allow_revcomp: bool,
    expected_strand: int,
    expected_start: int,
    expected_end: int,
):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _consensus("Alu", "SINE", "Alu")

    # Feature refers to AGP component_id "comp1"
    doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_feature(seq_region="comp1", s=10, e=20, strand=1, consensus=rc)],
    }

    json_path = _write_json(tmp_path / "in.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(json_path)])

    # Map comp1(1..100) onto chr1(100..199) with given orientation.
    agp = tmp_path / "test.agp"
    agp.write_text(
        f"chr1\t100\t199\t1\tW\tcomp1\t1\t100\t{orientation}\n",
        encoding="utf-8",
    )

    out_json = tmp_path / "out.json"

    combine_json.combine_repeat_json(
        json_manifest=manifest,
        out_json=out_json,
        chunk_re=CHUNK_RE,
        agp_file=agp,
        allow_revcomp=allow_revcomp,
        validate_inputs=False,
        validate_output=False,
    )

    out = json.loads(out_json.read_text(encoding="utf-8"))
    feat = out["repeat_features"][0]

    assert feat["seq_region"] == "chr1"
    assert feat["seq_region_start"] == expected_start
    assert feat["seq_region_end"] == expected_end
    assert feat["seq_region_strand"] == expected_strand


def test_combine_repeat_json_agp_driven_missing_component_raises(tmp_path: Path):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _consensus("Alu", "SINE", "Alu")

    doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_feature(seq_region="comp_missing", s=10, e=20, strand=1, consensus=rc)],
    }

    json_path = _write_json(tmp_path / "in.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(json_path)])

    agp = tmp_path / "test.agp"
    agp.write_text("chr1\t100\t199\t1\tW\tcomp1\t1\t100\t+\n", encoding="utf-8")

    with pytest.raises(KeyError, match=r"not found as an AGP component_id"):
        combine_json.combine_repeat_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
            agp_file=agp,
            allow_revcomp=False,
            validate_inputs=False,
            validate_output=False,
        )


def test_combine_repeat_json_missing_consensus_reference_raises(tmp_path: Path):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _consensus("X", "Y", "Z")

    doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [],
        "repeat_features": [_feature(seq_region="chr1", s=1, e=2, consensus=rc)],
    }

    p = _write_json(tmp_path / "a.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p)])

    with pytest.raises(ValueError, match=r"not present in repeat_consensus"):
        combine_json.combine_repeat_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
            validate_inputs=False,
            validate_output=False,
        )


def test_combine_repeat_json_analysis_mismatch_raises(tmp_path: Path):
    d1 = {
        "analysis": _analysis("a"),
        "source": _source("prov"),
        "repeat_features": [],
        "repeat_consensus": [],
    }
    d2 = {
        "analysis": _analysis("b"),
        "source": _source("prov"),
        "repeat_features": [],
        "repeat_consensus": [],
    }

    p1 = _write_json(tmp_path / "a.json", d1)
    p2 = _write_json(tmp_path / "b.json", d2)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p1), str(p2)])

    with pytest.raises(ValueError, match=r"Top-level 'analysis' differs"):
        combine_json.combine_repeat_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
            validate_inputs=False,
            validate_output=False,
        )


def test_combine_repeat_json_invalid_feature_strand_raises(tmp_path: Path):
    doc = {
        "analysis": _analysis("rm"),
        "source": _source("prov"),
        "repeat_consensus": [],
        "repeat_features": [
            {
                "seq_region": "chr1",
                "seq_region_start": 1,
                "seq_region_end": 2,
                "seq_region_strand": 0,  # invalid
                "repeat_start": 1,
                "repeat_end": 2,
            }
        ],
    }

    p = _write_json(tmp_path / "a.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p)])

    with pytest.raises(ValueError, match=r"repeat_feature\.seq_region_strand must be -1 or 1"):
        combine_json.combine_repeat_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
            validate_inputs=False,
            validate_output=False,
        )
