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

import argparse
import gzip
import hashlib
import json
import logging
from pathlib import Path
import re
from typing import cast

import pytest

from ensembl.io.genomio.features import combine_json
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
    return {
        "run_date": "2026-02-18T00:00:00Z",
        "logic_name": name,
        "display_label": name,
        "description": f"{name} analysis",
        "program": "test",
        "program_version": "0.0",
    }


def _source(provider: str = "prov") -> dict[str, object]:
    return {"source_provider": provider, "is_primary": True}


def _md5_key(rn: str, rc: str, rt: str, seq: str | None = None) -> str:
    norm = "".join((seq or "").split()).upper()
    payload = "\t".join([rn, rc, rt, norm])
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def _repeat_consensus(
    rn: str = "Alu",
    rc_class: str = "SINE",
    rt: str = "Alu",
    seq: str | None = None,
) -> combine_json.RepeatConsensus:
    rc: combine_json.RepeatConsensus = {
        "repeat_consensus_key": _md5_key(rn, rc_class, rt, seq),
        "repeat_name": rn,
        "repeat_class": rc_class,
        "repeat_type": rt,
    }
    if seq is not None:
        rc["repeat_consensus"] = seq
    return rc


def _repeat_feature(
    seq_region: str,
    s: int,
    e: int,
    strand: int = 1,
    consensus_key: str | None = None,
) -> combine_json.RepeatFeature:
    feat: combine_json.RepeatFeature = {
        "seq_region": seq_region,
        "seq_region_start": s,
        "seq_region_end": e,
        "seq_region_strand": strand,
        "repeat_start": 1,
        "repeat_end": 10,
    }
    if consensus_key is not None:
        feat["repeat_consensus"] = consensus_key
    return feat


def _ncrna_feature(seq_region: str, s: int, e: int, strand: int = 1) -> dict[str, object]:
    return {
        "seq_region": seq_region,
        "seq_region_start": s,
        "seq_region_end": e,
        "seq_region_strand": strand,
        "biotype": "miRNA",
        "score": 1.0,
        "target_name": "MIRTEST",
        "is_significant": True,
    }


def _trnascan_feature(seq_region: str, s: int, e: int, strand: int = 1) -> dict[str, object]:
    return {
        "seq_region": seq_region,
        "seq_region_start": s,
        "seq_region_end": e,
        "seq_region_strand": strand,
        "biotype": "tRNA",
        "score": 50.0,
        "isotype": "Phe",
        "anticodon": "GAA",
    }


def test_top_level_accumulator_get_required_raises_when_missing():
    acc = combine_json._TopLevelAccumulator()
    with pytest.raises(ValueError, match=r"Missing required top-level 'analysis'"):
        acc.get_required("analysis")


@pytest.mark.parametrize("gz", [False, True])
def test_load_json_document_accepts_object(tmp_path: Path, gz: bool):
    doc = {
        "analysis": _analysis(),
        "source": _source(),
        "repeat_features": [_repeat_feature(seq_region="chr1", s=1, e=2)],
    }
    p = tmp_path / ("in.json.gz" if gz else "in.json")
    (_write_json_gz if gz else _write_json)(p, doc)

    got = combine_json._load_json_document(p)
    assert isinstance(got, dict)
    assert "analysis" in got


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


def test_coerce_repeat_consensus_rejects_key_mismatch(tmp_path: Path):
    rc = _repeat_consensus("Alu", "SINE", "Alu", seq="ACGT")
    # Break it: change sequence but keep key
    bad = dict(rc)
    bad["repeat_consensus"] = "TTTT"

    with pytest.raises(ValueError, match=r"repeat_consensus_key mismatch"):
        combine_json._coerce_repeat_consensus(
            cast(combine_json.JsonValue, bad), source_path=tmp_path / "a.json"
        )


def test_merge_repeat_consensus_dedupes_identical(tmp_path: Path):
    rc = _repeat_consensus("Alu", "SINE", "Alu", seq="ACGT")
    combined: dict[str, combine_json.RepeatConsensus] = {}

    incoming = [cast(combine_json.JsonValue, rc)]

    combine_json._merge_repeat_consensus(combined, incoming, source_path=tmp_path / "a.json")
    combine_json._merge_repeat_consensus(combined, incoming, source_path=tmp_path / "b.json")

    assert len(combined) == 1


def test_coerce_repeat_feature_rejects_start_gt_end(tmp_path: Path):
    feat = _repeat_feature(seq_region="chr1", s=10, e=5)
    with pytest.raises(ValueError, match=r"seq_region_start > seq_region_end"):
        combine_json._coerce_repeat_feature(
            cast(combine_json.JsonValue, feat),
            source_path=tmp_path / "a.json",
        )


def test_coerce_ncrna_feature_rejects_start_gt_end(tmp_path: Path):
    feat = _ncrna_feature("chr1", 10, 5, strand=1)
    with pytest.raises(ValueError, match=r"seq_region_start > seq_region_end"):
        combine_json._coerce_ncrna_feature(
            cast(combine_json.JsonValue, feat),
            source_path=tmp_path / "a.json",
        )


@pytest.mark.parametrize(
    "name,s,e,expected_region,expected_start,expected_end",
    [
        ("chr1_chunk_start_11", 1, 5, "chr1", 11, 15),
        ("chr1", 1, 5, "chr1", 1, 5),
    ],
)
def test_lift_repeat_feature_seq_region_driven(
    tmp_path: Path,
    name: str,
    s: int,
    e: int,
    expected_region: str,
    expected_start: int,
    expected_end: int,
):
    feat = _repeat_feature(seq_region=name, s=s, e=e)

    out = cast(
        combine_json.RepeatFeature,
        combine_json._lift_feature_coords(
            feat,
            chunk_re=CHUNK_RE,
            agp_by_component=None,
            allow_revcomp=False,
            source_path=tmp_path / "x.json",
        ),
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

    forward_strand_feature = _repeat_feature(seq_region="compP", s=10, e=20, strand=1)
    forward_strand_result = cast(
        combine_json.RepeatFeature,
        combine_json._lift_feature_coords(
            forward_strand_feature,
            chunk_re=CHUNK_RE,
            agp_by_component={"compP": [forward_strand_entry]},
            allow_revcomp=False,
            source_path=tmp_path / "x.json",
        ),
    )
    assert forward_strand_result["seq_region"] == "objP"
    assert forward_strand_result["seq_region_start"] == 109
    assert forward_strand_result["seq_region_end"] == 119
    assert forward_strand_result["seq_region_strand"] == 1

    reverse_strand_feature = _repeat_feature(seq_region="compM", s=10, e=20, strand=1)
    reverse_strand_result = cast(
        combine_json.RepeatFeature,
        combine_json._lift_feature_coords(
            reverse_strand_feature,
            chunk_re=CHUNK_RE,
            agp_by_component={"compM": [reverse_strand_entry]},
            allow_revcomp=True,
            source_path=tmp_path / "x.json",
        ),
    )
    assert reverse_strand_result["seq_region"] == "objM"
    assert reverse_strand_result["seq_region_strand"] == -1


def test_lift_feature_coords_rejects_end_greater_than_start(tmp_path: Path):
    feat = _repeat_feature(seq_region="chr1", s=10, e=5)
    with pytest.raises(ValueError, match=r"seq_region_start > seq_region_end"):
        combine_json._lift_feature_coords(
            feat,
            chunk_re=CHUNK_RE,
            agp_by_component=None,
            allow_revcomp=False,
            source_path=tmp_path / "x.json",
        )


def test_detect_load_type_rejects_unknown_type(tmp_path: Path):
    with pytest.raises(ValueError, match=r"Cannot auto-detect schema kind"):
        combine_json._detect_load_type(
            {"unexpected_key": "value"},
            tmp_path / "x.json",
        )


def test_feature_consensus_key_returns_none_when_no_consensus(tmp_path: Path):
    feat = _repeat_feature(seq_region="chr1", s=1, e=10, consensus_key=None)
    assert combine_json._feature_consensus_key(feat, tmp_path / "x.json") is None


def test_feature_consensus_key_rejects_invalid_consensus_key(tmp_path: Path):
    feat = _repeat_feature(seq_region="chr1", s=1, e=10, consensus_key="invalid_key")
    with pytest.raises(ValueError, match=r"repeat_consensus must be a 32-char hex md5"):
        combine_json._feature_consensus_key(feat, tmp_path / "x.json")


def test_iterate_validated_documents_validate_false_skips_schema_validation(
    schema_validator_calls, tmp_path: Path
):
    doc = {
        "analysis": _analysis(),
        "source": _source(),
        "repeat_features": [_repeat_feature(seq_region="chr1", s=1, e=2)],
    }
    p = _write_json(tmp_path / "a.json", doc)

    got = list(combine_json._iterate_validated_documents([p], validate=False))
    assert got and got[0][0] == p
    assert schema_validator_calls == []


def test_write_and_validate_writes_newline_and_calls_schema_validator(schema_validator_calls, tmp_path: Path):
    out_json = tmp_path / "out.json"
    combine_json._write_and_validate(out_json, {"hello": "world"})

    assert out_json.read_text(encoding="utf-8").endswith("\n")
    assert schema_validator_calls, "schema_validator should have been called for the output"


def test_combine_repeat_json_paths_empty_manifest_raises(tmp_path: Path):
    with pytest.raises(ValueError, match=r"No JSON files were read from the manifest"):
        combine_json._combine_repeat_json_paths(
            json_paths=[],
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
            agp_by_component=None,
            allow_revcomp=False,
        )


def test_combine_feature_json_seq_region_driven_merge_and_coord_liftover(
    schema_validator_calls, tmp_path: Path
):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _repeat_consensus("Alu", "SINE", "Alu", "ACGT")
    rc_key = rc["repeat_consensus_key"]

    d1 = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_repeat_feature(seq_region="chr1_chunk_start_1", s=1, e=3, consensus_key=rc_key)],
    }
    d2 = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_repeat_feature(seq_region="chr1_chunk_start_4", s=1, e=2, consensus_key=rc_key)],
    }

    p1 = _write_json(tmp_path / "a.json", d1)
    p2 = _write_json(tmp_path / "b.json", d2)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p1), str(p2)])
    out_json = tmp_path / "out.json"

    combine_json.combine_feature_json(
        json_manifest=manifest,
        out_json=out_json,
        chunk_re=CHUNK_RE,
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
def test_combine_feature_json_agp_driven_lifts_coordinates_and_strand(
    tmp_path: Path,
    orientation: str,
    allow_revcomp: bool,
    expected_strand: int,
    expected_start: int,
    expected_end: int,
):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _repeat_consensus("Alu", "SINE", "Alu", "ACGT")
    rc_key = rc["repeat_consensus_key"]

    # Feature refers to AGP component_id "comp1"
    doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_repeat_feature(seq_region="comp1", s=10, e=20, strand=1, consensus_key=rc_key)],
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

    combine_json.combine_feature_json(
        json_manifest=manifest,
        out_json=out_json,
        chunk_re=CHUNK_RE,
        agp_file=agp,
        allow_revcomp=allow_revcomp,
    )

    out = json.loads(out_json.read_text(encoding="utf-8"))
    feat = out["repeat_features"][0]

    assert feat["seq_region"] == "chr1"
    assert feat["seq_region_start"] == expected_start
    assert feat["seq_region_end"] == expected_end
    assert feat["seq_region_strand"] == expected_strand


def test_combine_feature_json_rejects_empty_manifest(tmp_path: Path):
    manifest = write_manifest(tmp_path / "manifest.txt", [])
    with pytest.raises(ValueError, match=r"empty manifest"):
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
        )


def test_combine_feature_json_agp_driven_missing_component_raises(tmp_path: Path):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _repeat_consensus("Alu", "SINE", "Alu", "ACGT")
    rc_key = rc["repeat_consensus_key"]

    doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [
            _repeat_feature(seq_region="comp_missing", s=10, e=20, strand=1, consensus_key=rc_key)
        ],
    }

    json_path = _write_json(tmp_path / "in.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(json_path)])

    agp = tmp_path / "test.agp"
    agp.write_text("chr1\t100\t199\t1\tW\tcomp1\t1\t100\t+\n", encoding="utf-8")

    with pytest.raises(KeyError, match=r"not found as an AGP component_id"):
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
            agp_file=agp,
            allow_revcomp=False,
        )


def test_combine_feature_json_agp_driven_ambiguous_mapping_raises(tmp_path: Path):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _repeat_consensus("Alu", "SINE", "Alu", "ACGT")
    rc_key = rc["repeat_consensus_key"]

    # This feature range will match *both* AGP entries below (overlapping component spans).
    doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [rc],
        "repeat_features": [_repeat_feature(seq_region="comp1", s=60, e=70, strand=1, consensus_key=rc_key)],
    }

    json_path = _write_json(tmp_path / "in.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(json_path)])

    agp = tmp_path / "test.agp"
    agp.write_text(
        "\n".join(
            [
                # Both have component_id comp1, overlapping comp spans that contain 60..70
                "chr1\t100\t199\t1\tW\tcomp1\t1\t100\t+",
                "chr1\t300\t399\t2\tW\tcomp1\t50\t150\t+",
                "",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"Ambiguous AGP mapping"):
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
            agp_file=agp,
            allow_revcomp=False,
        )


def test_combine_feature_json_missing_consensus_reference_raises(tmp_path: Path):
    analysis = _analysis("rm")
    source = _source("prov")
    rc = _repeat_consensus("X", "Y", "Z", seq="ACGT")
    rc_key = rc["repeat_consensus_key"]

    doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [],
        "repeat_features": [_repeat_feature(seq_region="chr1", s=1, e=2, consensus_key=rc_key)],
    }

    p = _write_json(tmp_path / "a.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p)])

    with pytest.raises(ValueError, match=r"not present in repeat_consensus"):
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
        )


def test_combine_feature_json_analysis_mismatch_raises(tmp_path: Path):
    d1 = {
        "analysis": _analysis("a"),
        "source": _source("prov"),
        "repeat_features": [_repeat_feature(seq_region="chr1", s=1, e=2)],
        "repeat_consensus": [],
    }
    d2 = {
        "analysis": _analysis("b"),
        "source": _source("prov"),
        "repeat_features": [_repeat_feature(seq_region="chr1", s=1, e=2)],
        "repeat_consensus": [],
    }

    p1 = _write_json(tmp_path / "a.json", d1)
    p2 = _write_json(tmp_path / "b.json", d2)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p1), str(p2)])

    with pytest.raises(ValueError, match=r"Top-level 'analysis' differs"):
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
        )


def test_combine_ncrna_json_header_driven_merge_and_coord_liftover(schema_validator_calls, tmp_path: Path):
    analysis = _analysis("cmscan")
    source = _source("prov")

    d1 = {
        "analysis": analysis,
        "source": source,
        "ncrna_tool": "cmscan",
        "ncrna_features": [_ncrna_feature("chr1_chunk_start_1", 1, 3, strand=1)],
    }
    d2 = {
        "analysis": analysis,
        "source": source,
        "ncrna_tool": "cmscan",
        "ncrna_features": [_ncrna_feature("chr1_chunk_start_4", 1, 2, strand=1)],
    }

    p1 = _write_json(tmp_path / "a.json", d1)
    p2 = _write_json(tmp_path / "b.json", d2)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p1), str(p2)])

    out_json = tmp_path / "out.json"
    combine_json.combine_feature_json(json_manifest=manifest, out_json=out_json, chunk_re=CHUNK_RE)

    out = json.loads(out_json.read_text(encoding="utf-8"))
    assert out["analysis"] == analysis
    assert out["source"] == source
    assert out["ncrna_tool"] == "cmscan"
    assert len(out["ncrna_features"]) == 2
    assert out["ncrna_features"][0]["seq_region"] == "chr1"
    assert out["ncrna_features"][1]["seq_region_start"] == 4

    assert schema_validator_calls


@pytest.mark.parametrize(
    "orientation,allow_revcomp,expected_strand,expected_start,expected_end",
    [
        ("+", False, 1, 109, 119),
        ("-", True, -1, 180, 190),
    ],
)
def test_combine_ncrna_json_agp_driven_lifts_coordinates_and_strand(
    tmp_path: Path,
    orientation: str,
    allow_revcomp: bool,
    expected_strand: int,
    expected_start: int,
    expected_end: int,
):
    analysis = _analysis("cmscan")
    source = _source("prov")

    doc = {
        "analysis": analysis,
        "source": source,
        "ncrna_tool": "cmscan",
        "ncrna_features": [_ncrna_feature("comp1", 10, 20, strand=1)],
    }

    json_path = _write_json(tmp_path / "in.json", doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(json_path)])

    agp = tmp_path / "test.agp"
    agp.write_text(f"chr1\t100\t199\t1\tW\tcomp1\t1\t100\t{orientation}\n", encoding="utf-8")

    out_json = tmp_path / "out.json"
    combine_json.combine_feature_json(
        json_manifest=manifest,
        out_json=out_json,
        chunk_re=CHUNK_RE,
        agp_file=agp,
        allow_revcomp=allow_revcomp,
    )

    out = json.loads(out_json.read_text(encoding="utf-8"))
    feat = out["ncrna_features"][0]
    assert feat["seq_region"] == "chr1"
    assert feat["seq_region_start"] == expected_start
    assert feat["seq_region_end"] == expected_end
    assert feat["seq_region_strand"] == expected_strand


def test_combine_ncrna_json_tool_mismatch_raises(tmp_path: Path):
    analysis = _analysis("ncrna")
    source = _source("prov")

    d1 = {
        "analysis": analysis,
        "source": source,
        "ncrna_tool": "cmscan",
        "ncrna_features": [_ncrna_feature("chr1", 1, 2)],
    }
    d2 = {
        "analysis": analysis,
        "source": source,
        "ncrna_tool": "trnascan",
        "ncrna_features": [_trnascan_feature("chr1", 3, 4)],
    }

    p1 = _write_json(tmp_path / "a.json", d1)
    p2 = _write_json(tmp_path / "b.json", d2)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p1), str(p2)])

    with pytest.raises(ValueError, match=r"ncrna_tool.*differs"):
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
        )


def test_combine_feature_json_mixed_schema_kinds_raises(tmp_path: Path):
    analysis = _analysis("mix")
    source = _source("prov")

    repeat_doc = {
        "analysis": analysis,
        "source": source,
        "repeat_consensus": [],
        "repeat_features": [_repeat_feature("chr1", 1, 2)],
    }
    ncrna_doc = {
        "analysis": analysis,
        "source": source,
        "ncrna_tool": "cmscan",
        "ncrna_features": [_ncrna_feature("chr1", 1, 2)],
    }

    p1 = _write_json(tmp_path / "repeat.json", repeat_doc)
    p2 = _write_json(tmp_path / "ncrna.json", ncrna_doc)
    manifest = write_manifest(tmp_path / "manifest.txt", [str(p1), str(p2)])

    with pytest.raises(ValueError, match=r"Mixed load types detected"):
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=tmp_path / "out.json",
            chunk_re=CHUNK_RE,
        )


def test_parse_args_minimal_required_calls_init_logging(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    manifest = write_manifest(tmp_path / "manifest.txt", ["a.json"])
    out_json = tmp_path / "out.json"

    called: dict[str, object] = {}

    def fake_init_logging(args: argparse.Namespace) -> None:
        called["args"] = args

    monkeypatch.setattr(combine_json, "init_logging_with_args", fake_init_logging)

    args = combine_json.parse_args(
        [
            "--json-manifest",
            str(manifest),
            "--out-json",
            str(out_json),
        ]
    )

    assert args.json_manifest == manifest
    assert args.out_json == out_json
    assert not hasattr(args, "agp_file")
    assert args.allow_revcomp is False
    assert isinstance(args.chunk_id_regex, str)
    assert re.compile(args.chunk_id_regex).match("chr1_chunk_start_0")

    assert called["args"] is args


def test_parse_args_allow_revcomp_flag(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    manifest = write_manifest(tmp_path / "manifest.txt", ["a.json"])
    out_json = tmp_path / "out.json"

    monkeypatch.setattr(combine_json, "init_logging_with_args", lambda args: None)

    args = combine_json.parse_args(
        [
            "--json-manifest",
            str(manifest),
            "--out-json",
            str(out_json),
            "--allow-revcomp",
        ]
    )

    assert args.allow_revcomp is True


def test_parse_args_agp_file_sets_attribute(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    manifest = write_manifest(tmp_path / "manifest.txt", ["a.json"])
    out_json = tmp_path / "out.json"

    agp = tmp_path / "in.agp"
    agp.write_text("obj\t1\t1\t1\tW\tp1\t1\t1\t+\n")

    monkeypatch.setattr(combine_json, "init_logging_with_args", lambda args: None)

    args = combine_json.parse_args(
        [
            "--json-manifest",
            str(manifest),
            "--out-json",
            str(out_json),
            "--agp-file",
            str(agp),
        ]
    )

    assert args.agp_file == agp
    assert args.allow_revcomp is False


def test_parse_args_custom_chunk_id_regex(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    manifest = write_manifest(tmp_path / "manifest.txt", ["a.json"])
    out_json = tmp_path / "out.json"

    monkeypatch.setattr(combine_json, "init_logging_with_args", lambda args: None)

    custom = r"^(?P<base>.+)\.chunk\.(?P<start>\d+)$"

    args = combine_json.parse_args(
        [
            "--json-manifest",
            str(manifest),
            "--out-json",
            str(out_json),
            "--chunk-id-regex",
            custom,
        ]
    )

    assert args.chunk_id_regex == custom
    assert re.compile(args.chunk_id_regex).match("scaf_1.chunk.12")


def test_parse_args_missing_required_out_json_exits(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    manifest = write_manifest(tmp_path / "manifest.txt", ["a.json"])
    monkeypatch.setattr(combine_json, "init_logging_with_args", lambda args: None)

    with pytest.raises(SystemExit) as e:
        combine_json.parse_args(["--json-manifest", str(manifest)])

    assert e.value.code == 2


def test_parse_args_missing_manifest_file_exits(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    missing_manifest = tmp_path / "missing.txt"
    out_json = tmp_path / "out.json"
    monkeypatch.setattr(combine_json, "init_logging_with_args", lambda args: None)

    with pytest.raises(SystemExit) as e:
        combine_json.parse_args(["--json-manifest", str(missing_manifest), "--out-json", str(out_json)])

    assert e.value.code == 2


def test_parse_args_missing_agp_file_exits(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    manifest = write_manifest(tmp_path / "manifest.txt", ["a.json"])
    out_json = tmp_path / "out.json"
    missing_agp = tmp_path / "missing.agp"

    monkeypatch.setattr(combine_json, "init_logging_with_args", lambda args: None)

    with pytest.raises(SystemExit) as e:
        combine_json.parse_args(
            [
                "--json-manifest",
                str(manifest),
                "--out-json",
                str(out_json),
                "--agp-file",
                str(missing_agp),
            ]
        )

    assert e.value.code == 2


def test_main_logs_and_reraises_on_exception(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    manifest = write_manifest(tmp_path / "manifest.txt", ["a.json"])
    out_json = tmp_path / "out.json"

    monkeypatch.setattr(combine_json, "init_logging_with_args", lambda args: None)
    monkeypatch.setattr(combine_json, "validate_regex", lambda s: re.compile(s))

    def raise_main_exception(*args, **kwargs):
        raise RuntimeError("Simulated exception in main")

    monkeypatch.setattr(combine_json, "combine_feature_json", raise_main_exception)

    caplog.set_level(logging.ERROR)

    with pytest.raises(RuntimeError, match="Simulated exception in main"):
        combine_json.main(
            [
                "--json-manifest",
                str(manifest),
                "--out-json",
                str(out_json),
            ]
        )

    assert "Error combining feature JSON from files in" in caplog.text
    assert str(manifest) in caplog.text
