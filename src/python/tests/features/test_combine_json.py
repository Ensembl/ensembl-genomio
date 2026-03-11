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

from contextlib import nullcontext as does_not_raise
import hashlib
import json
from pathlib import Path
import re
from typing import cast, ContextManager
from unittest.mock import Mock, patch

from deepdiff import DeepDiff
import pytest
from pytest import param

from ensembl.io.genomio.features import combine_json


CHUNK_RE = re.compile(r"^(?P<base>.+)_chunk_start_(?P<start>\d+)$")


@pytest.fixture()
def schema_validator_calls(monkeypatch: pytest.MonkeyPatch):
    """
    Records calls made to schema_validator and replaces it with a no-op.

    Args:
        monkeypatch: Pytest monkeypatch fixture.

    Returns:
        List of ``(args, kwargs)`` tuples captured from calls to ``schema_validator``.
    """
    calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def _schema_validator_recorder(*args: object, **kwargs: object) -> None:
        calls.append((args, kwargs))

    monkeypatch.setattr(combine_json, "schema_validator", _schema_validator_recorder)
    return calls


def _analysis(name: str = "rm") -> dict[str, combine_json.JsonValue]:
    """
    Returns a standard analysis payload for assertions.

    Args:
        name: Analysis logic name.

    Returns:
        Analysis dictionary.
    """
    return {
        "run_date": "2026-02-18T00:00:00Z",
        "logic_name": name,
        "display_label": name,
        "description": f"{name} analysis",
        "program": "test",
        "program_version": "0.0",
    }


def _source(provider: str = "prov") -> dict[str, combine_json.JsonValue]:
    """
    Returns a standard source payload for assertions.

    Args:
        provider: Source provider string.

    Returns:
        Source dictionary.
    """
    return {"source_provider": provider, "is_primary": True}


def _md5_key(rn: str, rc_class: str, rt: str, seq: str | None = None) -> str:
    """
    Computes the repeat-consensus MD5 key used by the schema.

    Args:
        rn: Repeat name.
        rc_class: Repeat class.
        rt: Repeat type.
        seq: Optional consensus sequence.

    Returns:
        MD5 hash string.
    """
    norm = "".join((seq or "").split()).upper()
    payload = "\t".join([rn, rc_class, rt, norm])
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def _repeat_consensus(
    repeat_name: str = "Alu",
    repeat_class: str = "SINE",
    repeat_type: str = "Alu",
    sequence: str | None = None,
) -> combine_json.RepeatConsensus:
    """
    Builds a repeat-consensus object for tests.

    Args:
        repeat_name: Repeat name.
        repeat_class: Repeat class.
        repeat: Repeat type.
        sequence: Optional consensus sequence.

    Returns:
        Repeat-consensus dictionary.
    """
    rc: combine_json.RepeatConsensus = {
        "repeat_consensus_key": _md5_key(repeat_name, repeat_class, repeat_type, sequence),
        "repeat_name": repeat_name,
        "repeat_class": repeat_class,
        "repeat_type": repeat_type,
    }
    if sequence is not None:
        rc["repeat_consensus"] = sequence
    return rc


def _repeat_feature(
    seq_region: str,
    start: int,
    end: int,
    strand: int = 1,
    consensus_key: str | None = None,
) -> combine_json.RepeatFeature:
    """
    Builds a repeat-feature object for tests.

    Args:
        seq_region: Feature seq_region.
        start: Start coordinate.
        end: End coordinate.
        strand: Strand value.
        consensus_key: Optional repeat-consensus key.

    Returns:
        Repeat-feature dictionary.
    """
    feat: combine_json.RepeatFeature = {
        "seq_region": seq_region,
        "seq_region_start": start,
        "seq_region_end": end,
        "seq_region_strand": strand,
        "repeat_start": 1,
        "repeat_end": 10,
    }
    if consensus_key is not None:
        feat["repeat_consensus"] = consensus_key
    return feat


def _ncrna_feature(
    seq_region: str, start: int, end: int, strand: int = 1
) -> dict[str, combine_json.JsonValue]:
    """
    Builds a cmscan-style ncRNA feature object for tests.

    Args:
        seq_region: Feature seq_region.
        start: Start coordinate.
        end: End coordinate.
        strand: Strand value.

    Returns:
        ncRNA-feature dictionary.
    """
    return {
        "seq_region": seq_region,
        "seq_region_start": start,
        "seq_region_end": end,
        "seq_region_strand": strand,
        "biotype": "miRNA",
        "score": 1.0,
        "target_name": "MIRTEST",
        "is_significant": True,
    }


def _load_json(path: Path) -> dict[str, combine_json.JsonValue]:
    """
    Loads a JSON document from disk.

    Args:
        path: Path to JSON file.

    Returns:
        Parsed JSON document.
    """
    return cast(dict[str, combine_json.JsonValue], json.loads(path.read_text(encoding="utf-8")))


def test_top_level_accumulator_get_required_raises_when_missing():
    """Tests ``combine_json._TopLevelAccumulator.get_required()`` raises for a missing key."""
    acc = combine_json._TopLevelAccumulator()
    with pytest.raises(ValueError, match=r"Missing required top-level 'analysis'"):
        acc.get_required("analysis")


@pytest.mark.parametrize(
    "json_filename",
    [
        param("object.json", id="plain_json"),
        param("object.json.gz", id="gzipped_json"),
    ],
)
def test_load_json_document_accepts_object(data_dir: Path, json_filename: str) -> None:
    """
    Tests the `combine_json._load_json_document()` function for plain and gzipped JSON files.

    Args:
        data_dir: Module's test data directory fixture.
        json_filename: Name of input JSON file.
    """
    json_path = data_dir / "load_json" / json_filename

    document = combine_json._load_json_document(json_path)
    assert isinstance(document, dict)
    assert "analysis" in document
    assert "source" in document
    assert "repeat_features" in document


@pytest.mark.parametrize(
    "parts,start,end,expectation",
    [
        param(
            [
                combine_json.AgpEntry("obj", 1, 100, 1, "comp", 1, 100, "+"),
            ],
            10,
            20,
            does_not_raise(combine_json.AgpEntry("obj", 1, 100, 1, "comp", 1, 100, "+")),
            id="returns_matching_entry",
        ),
        param(
            [
                combine_json.AgpEntry("obj", 1, 100, 1, "comp", 1, 100, "+"),
                combine_json.AgpEntry("obj", 101, 200, 2, "comp", 101, 200, "+"),
            ],
            120,
            150,
            does_not_raise(combine_json.AgpEntry("obj", 101, 200, 2, "comp", 101, 200, "+")),
            id="returns_correct_entry_when_multiple_non_overlapping",
        ),
        param(
            [
                combine_json.AgpEntry("obj", 1, 100, 1, "comp", 1, 100, "+"),
            ],
            200,
            210,
            pytest.raises(ValueError, match=r"does not fit within any AGP span"),
            id="raises_when_no_span_matches",
        ),
        param(
            [
                combine_json.AgpEntry("obj", 1, 100, 1, "comp", 1, 100, "+"),
                combine_json.AgpEntry("obj", 1, 100, 2, "comp", 50, 150, "+"),
            ],
            60,
            70,
            pytest.raises(ValueError, match=r"Ambiguous AGP mapping"),
            id="raises_when_ambiguous",
        ),
    ],
)
def test_get_agp_entry_for_range(
    tmp_path: Path,
    parts: list[combine_json.AgpEntry],
    start: int,
    end: int,
    expectation: ContextManager,
) -> None:
    """
    Tests the `combine_json._get_agp_entry_for_range()` function.

    Args:
        tmp_path: Temporary directory provided by pytest.
        parts: AGP entries for a single component_id.
        start: Start coordinate on the component.
        end: End coordinate on the component.
        expectation: Context manager for the expected outcome of the test.
    """
    with expectation as expected:
        entry = combine_json._get_agp_entry_for_range(
            parts,
            start,
            end,
            component_id="comp",
            path=tmp_path / "x.agp",
        )
        assert entry == expected


@pytest.mark.parametrize(
    "value,coercer,expectation",
    [
        param(
            dict(_repeat_consensus("Alu", "SINE", "Alu", sequence="ACGT"), repeat_consensus="TTTT"),
            combine_json._coerce_repeat_consensus,
            pytest.raises(ValueError, match=r"repeat_consensus_key mismatch"),
            id="coerce_repeat_consensus_rejects_key_mismatch",
        ),
        param(
            _repeat_feature(seq_region="chr1", start=10, end=5),
            combine_json._coerce_repeat_feature,
            pytest.raises(ValueError, match=r"seq_region_start > seq_region_end"),
            id="coerce_repeat_feature_rejects_start_gt_end",
        ),
        param(
            _ncrna_feature("chr1", 10, 5, strand=1),
            combine_json._coerce_ncrna_feature,
            pytest.raises(ValueError, match=r"seq_region_start > seq_region_end"),
            id="coerce_ncrna_feature_rejects_start_gt_end",
        ),
    ],
)
def test_coercion_exceptions(
    tmp_path: Path,
    value: combine_json.JsonValue,
    coercer: combine_json.CoerceFeatureFn,
    expectation: ContextManager,
) -> None:
    """
    Tests exceptions raised by feature and consensus coercion helpers.

    Args:
        tmp_path: Temporary directory provided by pytest.
        value: Input object to coerce.
        coercer: Coercion function under test.
        expectation: Context manager for the expected exception.
    """
    with expectation:
        coercer(value, tmp_path / "a.json")


def test_merge_repeat_consensus_dedupes_identical(tmp_path: Path) -> None:
    """
    Tests ``combine_json._merge_repeat_consensus()`` deduplicates identical consensus entries.

    Args:
        tmp_path: Temporary directory provided by pytest.
    """
    rc = _repeat_consensus("Alu", "SINE", "Alu", sequence="ACGT")
    combined: dict[str, combine_json.RepeatConsensus] = {}
    incoming = [cast(combine_json.JsonValue, rc)]

    combine_json._merge_repeat_consensus(combined, incoming, source_path=tmp_path / "a.json")
    combine_json._merge_repeat_consensus(combined, incoming, source_path=tmp_path / "b.json")

    assert len(combined) == 1


@pytest.mark.parametrize(
    "feature,agp_by_component,allow_revcomp,expectation",
    [
        param(
            _repeat_feature(seq_region="chr1_chunk_start_11", start=1, end=5, strand=1),
            None,
            False,
            does_not_raise(
                {
                    "seq_region": "chr1",
                    "seq_region_start": 11,
                    "seq_region_end": 15,
                    "seq_region_strand": 1,
                }
            ),
            id="seq_region_driven_chunked",
        ),
        param(
            _repeat_feature(seq_region="chr1", start=1, end=5, strand=1),
            None,
            False,
            does_not_raise(
                {
                    "seq_region": "chr1",
                    "seq_region_start": 1,
                    "seq_region_end": 5,
                    "seq_region_strand": 1,
                }
            ),
            id="seq_region_driven_unchunked",
        ),
        param(
            _repeat_feature(seq_region="compP", start=10, end=20, strand=1),
            {
                "compP": [
                    combine_json.AgpEntry(
                        record="objP",
                        record_start=100,
                        record_end=199,
                        part_number=1,
                        part_id="compP",
                        part_start=1,
                        part_end=100,
                        orientation="+",
                    )
                ]
            },
            False,
            does_not_raise(
                {
                    "seq_region": "objP",
                    "seq_region_start": 109,
                    "seq_region_end": 119,
                    "seq_region_strand": 1,
                }
            ),
            id="agp_driven_forward",
        ),
        param(
            _repeat_feature(seq_region="compM", start=10, end=20, strand=1),
            {
                "compM": [
                    combine_json.AgpEntry(
                        record="objM",
                        record_start=100,
                        record_end=199,
                        part_number=1,
                        part_id="compM",
                        part_start=1,
                        part_end=100,
                        orientation="-",
                    )
                ]
            },
            True,
            does_not_raise(
                {
                    "seq_region": "objM",
                    "seq_region_start": 180,
                    "seq_region_end": 190,
                    "seq_region_strand": -1,
                }
            ),
            id="agp_driven_reverse",
        ),
        param(
            _repeat_feature(seq_region="chr1", start=10, end=5, strand=1),
            None,
            False,
            pytest.raises(ValueError, match=r"seq_region_start > seq_region_end"),
            id="rejects_end_less_than_start",
        ),
    ],
)
def test_lift_feature_coords(
    tmp_path: Path,
    feature: combine_json.RepeatFeature,
    agp_by_component: dict[str, list[combine_json.AgpEntry]] | None,
    allow_revcomp: bool,
    expectation: ContextManager,
) -> None:
    """
    Tests the `combine_json._lift_feature_coords()` function.

    Args:
        tmp_path: Temporary directory provided by pytest.
        feature: Input feature whose coordinates should be lifted.
        agp_by_component: Optional AGP component index used for AGP-driven lifting.
        allow_revcomp: Whether reverse-oriented AGP entries are permitted.
        expectation: Context manager for the expected outcome of the test.
    """
    with expectation as expected:
        out = cast(
            combine_json.RepeatFeature,
            combine_json._lift_feature_coords(
                feature,
                chunk_re=CHUNK_RE,
                agp_by_component=agp_by_component,
                allow_revcomp=allow_revcomp,
                source_path=tmp_path / "x.json",
            ),
        )
        assert out["seq_region"] == expected["seq_region"]
        assert out["seq_region_start"] == expected["seq_region_start"]
        assert out["seq_region_end"] == expected["seq_region_end"]
        assert out["seq_region_strand"] == expected["seq_region_strand"]


def test_detect_load_type_rejects_unknown_type(tmp_path: Path):
    """
    Tests ``combine_json._detect_load_type()`` rejects documents whose schema kind cannot be inferred.

    Args:
        tmp_path: Temporary directory provided by pytest.
    """
    with pytest.raises(ValueError, match=r"Cannot auto-detect schema kind"):
        combine_json._detect_load_type(
            {"unexpected_key": "value"},
            tmp_path / "x.json",
        )


@pytest.mark.parametrize(
    "feature,expectation",
    [
        param(
            _repeat_feature(seq_region="chr1", start=1, end=10, consensus_key=None),
            does_not_raise(None),
            id="returns_none_when_no_consensus",
        ),
        param(
            _repeat_feature(seq_region="chr1", start=1, end=10, consensus_key="A" * 32),
            does_not_raise("a" * 32),
            id="returns_lowercased_consensus_key",
        ),
        param(
            _repeat_feature(seq_region="chr1", start=1, end=10, consensus_key="invalid_key"),
            pytest.raises(ValueError, match=r"repeat_consensus must be a 32-char hex md5"),
            id="rejects_invalid_consensus_key",
        ),
    ],
)
def test_feature_consensus_key(
    tmp_path: Path,
    feature: combine_json.RepeatFeature,
    expectation: ContextManager,
) -> None:
    """
    Tests the `combine_json._feature_consensus_key()` function.

    Args:
        tmp_path: Temporary directory provided by pytest.
        feature: Repeat feature whose consensus key should be validated.
        expectation: Context manager for the expected outcome of the test.
    """
    with expectation as expected:
        assert combine_json._feature_consensus_key(feature, tmp_path / "x.json") == expected


def test_iterate_validated_documents_validate_false_skips_schema_validation(
    schema_validator_calls: list[tuple[tuple[object, ...], dict[str, object]]],
    data_dir: Path,
) -> None:
    """
    Tests ``combine_json._iterate_validated_documents()`` skips validation when requested.

    Args:
        schema_validator_calls: Captured ``schema_validator`` calls.
        data_dir: Module's test data directory fixture.
    """
    json_path = data_dir / "iterate_validated_documents" / "validate_false" / "a.json"

    documents = list(combine_json._iterate_validated_documents([json_path], validate=False))
    assert documents and documents[0][0] == json_path
    assert schema_validator_calls == []


def test_write_and_validate_writes_newline_and_calls_schema_validator(
    schema_validator_calls: list[tuple[tuple[object, ...], dict[str, object]]],
    tmp_path: Path,
) -> None:
    """
    Tests ``combine_json._write_and_validate()`` writes a trailing newline and validates output.

    Args:
        schema_validator_calls: Captured ``schema_validator`` calls.
        tmp_path: Temporary directory provided by pytest.
    """
    out_json = tmp_path / "out.json"
    combine_json._write_and_validate(out_json, {"hello": "world"})

    assert out_json.read_text(encoding="utf-8").endswith("\n")
    assert schema_validator_calls, "schema_validator should have been called for the output"


@pytest.mark.parametrize(
    "documents,feature_list_key,coerce_feature,agp_by_component,allow_revcomp,required_top_level_keys,expectation",
    [
        param(
            [
                (
                    Path("a.json"),
                    {
                        "analysis": _analysis("rm"),
                        "source": _source("prov"),
                        "repeat_features": [
                            _repeat_feature(seq_region="chr1_chunk_start_1", start=1, end=3, strand=1),
                        ],
                    },
                ),
                (
                    Path("b.json"),
                    {
                        "analysis": _analysis("rm"),
                        "source": _source("prov"),
                        "repeat_features": [
                            _repeat_feature(seq_region="chr1_chunk_start_4", start=1, end=2, strand=1),
                        ],
                    },
                ),
            ],
            "repeat_features",
            combine_json._coerce_repeat_feature,
            None,
            False,
            ["analysis", "source"],
            does_not_raise(
                {
                    "top_level": {
                        "analysis": _analysis("rm"),
                        "source": _source("prov"),
                    },
                    "features": [
                        _repeat_feature(seq_region="chr1", start=1, end=3, strand=1),
                        _repeat_feature(seq_region="chr1", start=4, end=5, strand=1),
                    ],
                    "nr_features": 2,
                }
            ),
            id="repeat_seq_region_driven_combines_and_lifts",
        ),
        param(
            [
                (
                    Path("a.json"),
                    {
                        "analysis": _analysis("cmscan"),
                        "source": _source("prov"),
                        "ncrna_tool": "cmscan",
                        "ncrna_features": [
                            _ncrna_feature(seq_region="chr1_chunk_start_1", start=1, end=3, strand=1),
                        ],
                    },
                ),
                (
                    Path("b.json"),
                    {
                        "analysis": _analysis("cmscan"),
                        "source": _source("prov"),
                        "ncrna_tool": "cmscan",
                        "ncrna_features": [
                            _ncrna_feature(seq_region="chr1_chunk_start_4", start=1, end=2, strand=1),
                        ],
                    },
                ),
            ],
            "ncrna_features",
            combine_json._coerce_ncrna_feature,
            None,
            False,
            ["analysis", "source", "ncrna_tool"],
            does_not_raise(
                {
                    "top_level": {
                        "analysis": _analysis("cmscan"),
                        "source": _source("prov"),
                        "ncrna_tool": "cmscan",
                    },
                    "features": [
                        _ncrna_feature(seq_region="chr1", start=1, end=3, strand=1),
                        _ncrna_feature(seq_region="chr1", start=4, end=5, strand=1),
                    ],
                    "nr_features": 2,
                }
            ),
            id="ncrna_seq_region_driven_combines_and_lifts",
        ),
        param(
            [
                (
                    Path("a.json"),
                    {
                        "analysis": _analysis("rm"),
                        "source": _source("prov"),
                        "repeat_features": [
                            _repeat_feature(seq_region="comp1", start=10, end=20, strand=1),
                        ],
                    },
                ),
            ],
            "repeat_features",
            combine_json._coerce_repeat_feature,
            {
                "comp1": [
                    combine_json.AgpEntry(
                        record="chr1",
                        record_start=100,
                        record_end=199,
                        part_number=1,
                        part_id="comp1",
                        part_start=1,
                        part_end=100,
                        orientation="+",
                    )
                ]
            },
            False,
            ["analysis", "source"],
            does_not_raise(
                {
                    "top_level": {
                        "analysis": _analysis("rm"),
                        "source": _source("prov"),
                    },
                    "features": [
                        _repeat_feature(seq_region="chr1", start=109, end=119, strand=1),
                    ],
                    "nr_features": 1,
                }
            ),
            id="repeat_agp_driven_forward_lifts",
        ),
        param(
            [
                (
                    Path("a.json"),
                    {
                        "analysis": _analysis("cmscan"),
                        "source": _source("prov"),
                        "ncrna_tool": "cmscan",
                        "ncrna_features": [
                            _ncrna_feature(seq_region="comp1", start=10, end=20, strand=1),
                        ],
                    },
                ),
            ],
            "ncrna_features",
            combine_json._coerce_ncrna_feature,
            {
                "comp1": [
                    combine_json.AgpEntry(
                        record="chr1",
                        record_start=100,
                        record_end=199,
                        part_number=1,
                        part_id="comp1",
                        part_start=1,
                        part_end=100,
                        orientation="-",
                    )
                ]
            },
            True,
            ["analysis", "source", "ncrna_tool"],
            does_not_raise(
                {
                    "top_level": {
                        "analysis": _analysis("cmscan"),
                        "source": _source("prov"),
                        "ncrna_tool": "cmscan",
                    },
                    "features": [
                        _ncrna_feature(seq_region="chr1", start=180, end=190, strand=-1),
                    ],
                    "nr_features": 1,
                }
            ),
            id="ncrna_agp_driven_reverse_lifts_and_flips_strand",
        ),
        param(
            [],
            "repeat_features",
            combine_json._coerce_repeat_feature,
            None,
            False,
            ["analysis", "source"],
            pytest.raises(ValueError, match=r"No JSON files were read from the manifest"),
            id="empty_documents_raises",
        ),
        param(
            [
                (
                    Path("a.json"),
                    {
                        "analysis": _analysis("a"),
                        "source": _source("prov"),
                        "repeat_features": [
                            _repeat_feature(seq_region="chr1", start=1, end=2, strand=1),
                        ],
                    },
                ),
                (
                    Path("b.json"),
                    {
                        "analysis": _analysis("b"),
                        "source": _source("prov"),
                        "repeat_features": [
                            _repeat_feature(seq_region="chr1", start=3, end=4, strand=1),
                        ],
                    },
                ),
            ],
            "repeat_features",
            combine_json._coerce_repeat_feature,
            None,
            False,
            ["analysis", "source"],
            pytest.raises(ValueError, match=r"Top-level 'analysis' differs"),
            id="analysis_mismatch_raises",
        ),
        param(
            [
                (
                    Path("a.json"),
                    {
                        "analysis": _analysis("cmscan"),
                        "source": _source("prov"),
                        "ncrna_tool": "cmscan",
                        "ncrna_features": [
                            _ncrna_feature(seq_region="chr1", start=1, end=2, strand=1),
                        ],
                    },
                ),
                (
                    Path("b.json"),
                    {
                        "analysis": _analysis("cmscan"),
                        "source": _source("prov"),
                        "ncrna_tool": "trnascan",
                        "ncrna_features": [
                            _ncrna_feature(seq_region="chr1", start=3, end=4, strand=1),
                        ],
                    },
                ),
            ],
            "ncrna_features",
            combine_json._coerce_ncrna_feature,
            None,
            False,
            ["analysis", "source", "ncrna_tool"],
            pytest.raises(ValueError, match=r"Top-level 'ncrna_tool' differs"),
            id="ncrna_tool_mismatch_raises",
        ),
    ],
)
def test_combine_feature_docs(
    tmp_path: Path,
    documents: list[tuple[Path, dict[str, combine_json.JsonValue]]],
    feature_list_key: str,
    coerce_feature: combine_json.CoerceFeatureFn,
    agp_by_component: dict[str, list[combine_json.AgpEntry]] | None,
    allow_revcomp: bool,
    required_top_level_keys: list[str],
    expectation: ContextManager,
) -> None:
    """
    Tests the `combine_json._combine_feature_docs()` function.

    Args:
        tmp_path: Temporary directory provided by pytest.
        documents: Input `(path, document)` pairs.
        feature_list_key: Top-level key containing the feature list.
        coerce_feature: Feature coercion function under test.
        agp_by_component: Optional AGP component index used for AGP-driven lifting.
        allow_revcomp: Boolean indicating whether minus strand entries are permitted.
        required_top_level_keys: Top-level keys that must be identical across all inputs.
        expectation: Context manager for the expected outcome of the test.
    """
    resolved_documents = [(tmp_path / path.name, document) for path, document in documents]

    with expectation as expected:
        top_level, features, nr_features = combine_json._combine_feature_docs(
            documents=resolved_documents,
            feature_list_key=feature_list_key,
            coerce_feature=coerce_feature,
            chunk_re=CHUNK_RE,
            agp_by_component=agp_by_component,
            allow_revcomp=allow_revcomp,
            required_top_level_keys=required_top_level_keys,
        )
        assert top_level == expected["top_level"]
        assert features == expected["features"]
        assert nr_features == expected["nr_features"]


@pytest.mark.parametrize(
    "test_dir_name,agp_filename,allow_revcomp,expectation",
    [
        param(
            "seq_region",
            None,
            False,
            does_not_raise(
                {
                    "analysis": _analysis("rm"),
                    "source": _source("prov"),
                    "repeat_consensus_len": 1,
                    "repeat_features": [
                        _repeat_feature(
                            seq_region="chr1",
                            start=1,
                            end=3,
                            strand=1,
                            consensus_key=_md5_key("Alu", "SINE", "Alu", "ACGT"),
                        ),
                        _repeat_feature(
                            seq_region="chr1",
                            start=4,
                            end=5,
                            strand=1,
                            consensus_key=_md5_key("Alu", "SINE", "Alu", "ACGT"),
                        ),
                    ],
                }
            ),
            id="seq_region_driven_repeat_merge_and_liftover",
        ),
        param(
            "agp_forward",
            "test.agp",
            False,
            does_not_raise(
                {
                    "analysis": _analysis("rm"),
                    "source": _source("prov"),
                    "repeat_consensus_len": 1,
                    "repeat_features": [
                        _repeat_feature(
                            seq_region="chr1",
                            start=109,
                            end=119,
                            strand=1,
                            consensus_key=_md5_key("Alu", "SINE", "Alu", "ACGT"),
                        ),
                    ],
                }
            ),
            id="agp_driven_forward_repeat_merge_and_liftover",
        ),
        param(
            "agp_reverse",
            "test.agp",
            True,
            does_not_raise(
                {
                    "analysis": _analysis("rm"),
                    "source": _source("prov"),
                    "repeat_consensus_len": 1,
                    "repeat_features": [
                        _repeat_feature(
                            seq_region="chr1",
                            start=180,
                            end=190,
                            strand=-1,
                            consensus_key=_md5_key("Alu", "SINE", "Alu", "ACGT"),
                        ),
                    ],
                }
            ),
            id="agp_driven_reverse_repeat_merge_and_liftover",
        ),
        param(
            "missing_consensus",
            None,
            False,
            pytest.raises(ValueError, match=r"not present in repeat_consensus"),
            id="missing_consensus_reference_raises",
        ),
    ],
)
def test_combine_repeat_json_paths(
    data_dir: Path,
    tmp_path: Path,
    test_dir_name: str,
    agp_filename: str | None,
    allow_revcomp: bool,
    expectation: ContextManager,
) -> None:
    """
    Tests the `combine_json._combine_repeat_json_paths()` function.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        test_dir_name: Name of data subdirectory for the test case.
        agp_filename: Optional AGP filename.
        allow_revcomp: Boolean indicating whether minus strand entries are permitted.
        expectation: Context manager for the expected outcome of the test.
    """
    test_dir = data_dir / "combine_repeat" / test_dir_name
    json_paths = combine_json.get_paths_from_manifest(test_dir / "manifest.txt")
    agp_file = test_dir / agp_filename if agp_filename is not None else None
    agp_by_component = None
    if agp_file is not None:
        agp_by_component = combine_json.build_component_index(combine_json.parse_agp(agp_file, allow_revcomp))

    out_json = tmp_path / f"{test_dir_name}.out.json"

    with expectation as expected:
        combine_json._combine_repeat_json_paths(
            json_paths=json_paths,
            out_json=out_json,
            chunk_re=CHUNK_RE,
            agp_by_component=agp_by_component,
            allow_revcomp=allow_revcomp,
        )

        out = _load_json(out_json)
        assert out["analysis"] == expected["analysis"]
        assert out["source"] == expected["source"]
        assert len(cast(list[object], out["repeat_consensus"])) == expected["repeat_consensus_len"]
        assert out["repeat_features"] == expected["repeat_features"]


@pytest.mark.parametrize(
    "test_dir_name,agp_filename,allow_revcomp,expectation",
    [
        param(
            "seq_region",
            None,
            False,
            does_not_raise(
                {
                    "analysis": _analysis("cmscan"),
                    "source": _source("prov"),
                    "ncrna_tool": "cmscan",
                    "ncrna_features": [
                        _ncrna_feature("chr1", 1, 3, strand=1),
                        _ncrna_feature("chr1", 4, 5, strand=1),
                        _ncrna_feature("chr1", 6, 7, strand=1),
                    ],
                }
            ),
            id="seq_region_driven_ncrna_merge_and_liftover",
        ),
        param(
            "agp_forward",
            "test.agp",
            False,
            does_not_raise(
                {
                    "analysis": _analysis("cmscan"),
                    "source": _source("prov"),
                    "ncrna_tool": "cmscan",
                    "ncrna_features": [
                        _ncrna_feature("chr1", 109, 119, strand=1),
                    ],
                }
            ),
            id="agp_driven_forward_ncrna_merge_and_liftover",
        ),
        param(
            "agp_reverse",
            "test.agp",
            True,
            does_not_raise(
                {
                    "analysis": _analysis("cmscan"),
                    "source": _source("prov"),
                    "ncrna_tool": "cmscan",
                    "ncrna_features": [
                        _ncrna_feature("chr1", 180, 190, strand=-1),
                    ],
                }
            ),
            id="agp_driven_reverse_ncrna_merge_and_liftover",
        ),
    ],
)
def test_combine_ncrna_json_paths(
    data_dir: Path,
    tmp_path: Path,
    test_dir_name: str,
    agp_filename: str | None,
    allow_revcomp: bool,
    expectation: ContextManager,
) -> None:
    """
    Tests the `combine_json._combine_ncrna_json_paths()` function.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        test_dir_name: Name of data subdirectory for the test case.
        agp_filename: Optional AGP filename.
        allow_revcomp: Boolean indicating whether minus strand entries are permitted.
        expectation: Context manager for the expected outcome of the test.
    """
    test_dir = data_dir / "combine_ncrna" / test_dir_name
    json_paths = combine_json.get_paths_from_manifest(test_dir / "manifest.txt")
    agp_file = test_dir / agp_filename if agp_filename is not None else None
    agp_by_component = None
    if agp_file is not None:
        agp_by_component = combine_json.build_component_index(combine_json.parse_agp(agp_file, allow_revcomp))

    out_json = tmp_path / f"{test_dir_name}.out.json"

    with expectation as expected:
        combine_json._combine_ncrna_json_paths(
            json_paths=json_paths,
            out_json=out_json,
            chunk_re=CHUNK_RE,
            agp_by_component=agp_by_component,
            allow_revcomp=allow_revcomp,
        )

        out = _load_json(out_json)
        assert out["analysis"] == expected["analysis"]
        assert out["source"] == expected["source"]
        assert out["ncrna_tool"] == expected["ncrna_tool"]
        assert out["ncrna_features"] == expected["ncrna_features"]


@pytest.mark.parametrize(
    "test_dir_name,agp_filename,allow_revcomp,expectation",
    [
        param(
            "combine_ncrna/seq_region",
            None,
            False,
            does_not_raise(
                {
                    "analysis": _analysis("cmscan"),
                    "source": _source("prov"),
                    "ncrna_tool": "cmscan",
                    "ncrna_features": [
                        _ncrna_feature("chr1", 1, 3, strand=1),
                        _ncrna_feature("chr1", 4, 5, strand=1),
                        _ncrna_feature("chr1", 6, 7, strand=1),
                    ],
                }
            ),
            id="ncrna_success",
        ),
        param(
            "empty_manifest",
            None,
            False,
            pytest.raises(ValueError, match=r"empty manifest"),
            id="rejects_empty_manifest",
        ),
        param(
            "combine_repeat/missing_component",
            "test.agp",
            False,
            pytest.raises(KeyError, match=r"not found as an AGP component_id"),
            id="agp_driven_missing_component_raises",
        ),
        param(
            "combine_repeat/ambiguous_mapping",
            "test.agp",
            False,
            pytest.raises(ValueError, match=r"Ambiguous AGP mapping"),
            id="agp_driven_ambiguous_mapping_raises",
        ),
        param(
            "mixed_schema_kinds",
            None,
            False,
            pytest.raises(ValueError, match=r"Mixed load types detected"),
            id="mixed_schema_kinds_raises",
        ),
    ],
)
def test_combine_feature_json(
    data_dir: Path,
    tmp_path: Path,
    test_dir_name: str,
    agp_filename: str | None,
    allow_revcomp: bool,
    expectation: ContextManager,
) -> None:
    """
    Tests ``combine_json.combine_feature_json()`` for successful and failing file-based inputs.

    Args:
        data_dir: Module's test data directory fixture.
        tmp_path: Temporary directory provided by pytest.
        test_dir_name: Name of data subdirectory for the test case.
        agp_filename: Optional AGP filename.
        allow_revcomp: Boolean indicating whether minus strand entries are permitted.
        expectation: Context manager for the expected outcome of the test.
    """
    test_dir = data_dir / test_dir_name
    manifest = test_dir / "manifest.txt"
    agp_file = test_dir / agp_filename if agp_filename is not None else None
    out_json = tmp_path / "out.json"

    with expectation as expected:
        combine_json.combine_feature_json(
            json_manifest=manifest,
            out_json=out_json,
            chunk_re=CHUNK_RE,
            agp_file=agp_file,
            allow_revcomp=allow_revcomp,
        )

        out = _load_json(out_json)
        assert out["analysis"] == expected["analysis"]
        assert out["source"] == expected["source"]
        assert out["ncrna_tool"] == expected["ncrna_tool"]
        assert out["ncrna_features"] == expected["ncrna_features"]


@pytest.mark.parametrize(
    "arg_list, expected_params",
    [
        param(
            ["--json-manifest", __file__, "--out-json", "out.json"],
            {
                "json_manifest": str(__file__),
                "out_json": "out.json",
                "allow_revcomp": False,
                "chunk_id_regex": CHUNK_RE.pattern,
                "log_level": "WARNING",
            },
            id="Default args",
        ),
        param(
            [
                "--json-manifest",
                __file__,
                "--out-json",
                "out.json",
                "--agp-file",
                __file__,
                "--allow-revcomp",
                "--chunk-id-regex",
                r"^(?P<base>.+)_start_(?P<start>\d+)$",
            ],
            {
                "json_manifest": str(__file__),
                "out_json": "out.json",
                "agp_file": str(__file__),
                "allow_revcomp": True,
                "chunk_id_regex": r"^(?P<base>.+)_start_(?P<start>\d+)$",
                "log_level": "WARNING",
            },
            id="New arg values",
        ),
    ],
)
def test_parse_args(arg_list: list[str], expected_params: dict[str, str | bool]) -> None:
    """
    Tests the `combine_json.parse_args()` function.

    Args:
        arg_list: List of command line arguments to parse.
        expected_params: Expected dictionary of parameters from parsing arguments.
    """
    args = combine_json.parse_args(arg_list)

    # DeepDiff cannot compare Path objects
    setattr(args, "json_manifest", str(args.json_manifest))
    setattr(args, "out_json", str(args.out_json))
    if hasattr(args, "agp_file"):
        setattr(args, "agp_file", str(args.agp_file))

    assert not DeepDiff(vars(args), expected_params)


@patch("ensembl.io.genomio.features.combine_json.combine_feature_json")
def test_main(mock_combine_feature_json: Mock, tmp_path: Path) -> None:
    """
    Tests the `combine_json.main()` function (entry point).

    Args:
        mock_combine_feature_json: Mock object for the `combine_json.combine_feature_json()` function.
        tmp_path: Temporary directory provided by pytest.
    """
    manifest_path = tmp_path / "manifest.txt"
    out_json = tmp_path / "out.json"
    manifest_path.touch()

    combine_json.main(["--json-manifest", str(manifest_path), "--out-json", str(out_json)])

    mock_combine_feature_json.assert_called_once_with(
        json_manifest=manifest_path,
        out_json=out_json,
        chunk_re=CHUNK_RE,
        agp_file=None,
        allow_revcomp=False,
    )


@patch("ensembl.io.genomio.features.combine_json.combine_feature_json")
def test_main_raises_exception(mock_combine_feature_json: Mock, tmp_path: Path) -> None:
    """
    Tests the `combine_json.main()` function (entry point) raises exception on error.

    Args:
        mock_combine_feature_json: Mock object for the `combine_json.combine_feature_json()` function.
        tmp_path: Temporary directory provided by pytest.
    """
    manifest_path = tmp_path / "manifest.txt"
    out_json = tmp_path / "out.json"
    manifest_path.touch()

    mock_combine_feature_json.side_effect = RuntimeError("Mocked exception")

    with pytest.raises(RuntimeError, match=r"Mocked exception"):
        combine_json.main(["--json-manifest", str(manifest_path), "--out-json", str(out_json)])
