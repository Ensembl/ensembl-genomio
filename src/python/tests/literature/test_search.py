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
"""Unit testing of `ensembl.io.genomio.literature.search` module.

Typical usage example::
    $ pytest test_search.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from rank_bm25 import BM25Okapi

from ensembl.io.genomio.literature import search


@pytest.mark.parametrize(
    ("section", "expected"),
    [
        pytest.param("Abstract", 2.0, id="Abstract"),
        pytest.param("2. Materials and Methods", 0.8, id="Methods"),
        pytest.param("Author Contributions", 0.0, id="Noise section"),
        pytest.param("Table 1", 1.0, id="Unmatched"),
    ],
)
def test_get_section_weight(section: str, expected: float) -> None:
    """Test :func:`search.get_section_weight`.

    Args:
        section: Section title.
        expected: Expected weight.

    """
    assert search.get_section_weight(section) == expected


@patch("ensembl.io.genomio.literature.search._embedder", None)
@patch("ensembl.io.genomio.literature.search.SentenceTransformer")
def test_get_embedder(mock_model: MagicMock) -> None:
    """Test :func:`search.get_embedder` loads the model once and reuses it.

    Args:
        mock_model: Patched `SentenceTransformer` class.

    """
    assert search.get_embedder() is mock_model.return_value
    assert search.get_embedder() is mock_model.return_value
    mock_model.assert_called_once()


@patch("ensembl.io.genomio.literature.search.get_embedder")
def test_build_index(mock_embedder: MagicMock, tmp_path: Path) -> None:
    """Test :func:`search.build_index` encodes, tokenises and persists the chunks.

    Args:
        mock_embedder: Patched `get_embedder()`.
        tmp_path: Pytest-provided temporary directory.

    """
    mock_embedder.return_value.encode.return_value = [[1.0, 0.0], [0.0, 1.0]]
    chunks = [
        {"section": "Abstract", "text": "Peanut is Tetraploid"},
        {"section": "Methods", "text": "2n = 40"},
    ]
    index_dir = tmp_path / "index"
    with patch.object(search, "INDEX_DIR", str(index_dir)):
        embeddings, bm25, meta, texts = search.build_index(chunks, paper_id="PMC1")
    assert mock_embedder.return_value.encode.call_args.args[0] == [
        "passage: Peanut is Tetraploid",
        "passage: 2n = 40",
    ]
    assert embeddings.tolist() == [[1.0, 0.0], [0.0, 1.0]]
    assert isinstance(bm25, BM25Okapi)
    assert (meta["paper_id"].tolist(), texts) == (["PMC1", "PMC1"], ["Peanut is Tetraploid", "2n = 40"])
    assert json.loads((index_dir / "bm25_tokens.json").read_text()) == [
        ["peanut", "is", "tetraploid"],
        ["2n", "=", "40"],
    ]


def test_load_index(tmp_path: Path) -> None:
    """Test :func:`search.load_index` reads back a persisted index.

    Args:
        tmp_path: Pytest-provided temporary directory.

    """
    np.save(tmp_path / "embeddings.npy", np.array([[1.0, 0.0]], dtype="float32"))
    (tmp_path / "bm25_tokens.json").write_text(json.dumps([["peanut", "tetraploid"]]))
    (tmp_path / "metadata.json").write_text(
        json.dumps([{"section": "Abstract", "text": "Peanut tetraploid"}])
    )
    with patch.object(search, "INDEX_DIR", str(tmp_path)):
        embeddings, bm25, meta, texts = search.load_index()
    assert embeddings.tolist() == [[1.0, 0.0]]
    assert isinstance(bm25, BM25Okapi)
    assert (meta["section"].tolist(), texts) == (["Abstract"], ["Peanut tetraploid"])


def test_rrf_fuse() -> None:
    """Test :func:`search.rrf_fuse`."""
    fused = search.rrf_fuse(np.array([3.0, 1.0, 2.0]), np.array([0.1, 0.5, 0.9]), k=2, rrf_k=1)
    assert [index for index, _ in fused] == [2, 0]
    assert [score for _, score in fused] == pytest.approx([1 / 2 + 1 / 1, 1 / 1 + 1 / 3])


@pytest.mark.parametrize(
    ("top_k", "expected_sections"),
    [
        pytest.param(2, ["Results", "Abstract"], id="Re-ranked and truncated"),
        pytest.param(5, ["Results", "Abstract", "Methods"], id="Noise section dropped"),
    ],
)
@patch("ensembl.io.genomio.literature.search.rrf_fuse")
@patch("ensembl.io.genomio.literature.search.get_embedder")
def test_search(
    mock_embedder: MagicMock, mock_fuse: MagicMock, top_k: int, expected_sections: list[str]
) -> None:
    """Test :func:`search.search` drops noise sections, weights and re-ranks the fused results.

    Args:
        mock_embedder: Patched `get_embedder()`.
        mock_fuse: Patched `rrf_fuse()`.
        top_k: Number of results to return.
        expected_sections: Sections of the expected results, in rank order.

    """
    mock_embedder.return_value.encode.return_value = np.array([[1.0, 0.0]])
    bm25 = MagicMock(**{"get_scores.return_value": [0.1, 0.2, 0.3, 0.4]})
    embeddings = np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.5], [0.2, 0.8]])
    meta = pd.DataFrame(
        {"section": ["Abstract", "Data availability", "Methods", "Results"], "paper_id": ["P"] * 4}
    )
    texts = ["t0", "t1", "t2", "t3"]
    # Fused order: Methods, Abstract, noise, Results; weighting moves Results to the top.
    mock_fuse.return_value = [(2, 1.0), (0, 0.5), (1, 0.9), (3, 0.9)]

    results = search.search("Ploidy Level", embeddings, bm25, meta, texts, top_k=top_k)

    bm25.get_scores.assert_called_once_with(["ploidy", "level"])
    assert mock_embedder.return_value.encode.call_args.args[0] == ["query: Ploidy Level"]
    assert mock_fuse.call_args.args[1].tolist() == [1.0, 0.0, 0.5, 0.2]
    all_results = {
        "Results": {"score": 1.08, "raw_score": 0.9, "section": "Results", "text": "t3", "paper_id": "P"},
        "Abstract": {"score": 1.0, "raw_score": 0.5, "section": "Abstract", "text": "t0", "paper_id": "P"},
        "Methods": {"score": 0.8, "raw_score": 1.0, "section": "Methods", "text": "t2", "paper_id": "P"},
    }
    assert results == [{"rank": i, **all_results[section]} for i, section in enumerate(expected_sections, 1)]


@patch("ensembl.io.genomio.literature.search.resolve_ploidy_fields")
def test_infer_ploidy(mock_resolve: MagicMock) -> None:
    """Test :func:`search.infer_ploidy` turns weighted results into resolver segments.

    Args:
        mock_resolve: Patched `resolve_ploidy_fields()`.

    """
    results = [
        {"section": "Abstract", "text": "tetraploid"},
        {"section": "Author contributions", "text": "diploid"},
        {"section": "Results", "text": "2n = 4x = 40"},
    ]
    assert search.infer_ploidy(results, chromosome_count=20, names=["arachis"]) is mock_resolve.return_value
    mock_resolve.assert_called_once_with(
        [("tetraploid", 2.0, True), ("2n = 4x = 40", 1.2, False)], chromosome_count=20, names=["arachis"]
    )


@pytest.mark.parametrize(
    ("results", "expected"),
    [
        pytest.param(
            [
                {"section": "Data availability", "text": "female donor"},
                {"section": "Methods", "text": "A single female plant was sequenced."},
            ],
            {"value": "female", "section": "Methods", "context": "A single female plant was sequenced."},
            id="Noise section skipped",
        ),
        pytest.param(
            [{"section": "Results", "text": "The ZW system"}],
            {"value": "ZW", "section": "Results", "context": "The ZW system"},
            id="Sex chromosome system",
        ),
        pytest.param(
            [{"section": "Results", "text": "Malesherbia"}],
            {"value": "unknown", "evidence": None},
            id="No keyword",
        ),
    ],
)
def test_infer_sex(results: list[dict], expected: dict) -> None:
    """Test :func:`search.infer_sex`.

    Args:
        results: Ranked search results.
        expected: Expected sex call.

    """
    assert search.infer_sex(results) == expected


@pytest.mark.xfail(strict=True, reason="search.py:248 finds the context with str.find(), not the regex match")
def test_infer_sex_context_surrounds_match() -> None:
    """Test :func:`search.infer_sex` returns context around the matched word, not an earlier substring."""
    text = (
        "Leaves from several females were pooled for RNA. "
        + "x" * 80
        + " The genome of a single male was sequenced."
    )
    assert "single male" in search.infer_sex([{"section": "Methods", "text": text}])["context"]


@pytest.mark.parametrize(
    ("results", "expected"),
    [
        pytest.param(
            [
                {"section": "Results", "text": "a", "score": 2.0},
                {"section": "Introduction", "text": "b", "score": 1.0},
            ],
            {"text": "b", "section": "Introduction", "score": 1.0},
            id="Priority section preferred",
        ),
        pytest.param(
            [{"section": "Results", "text": "a", "score": 2.0}],
            {"text": "a", "section": "Results", "score": 2.0},
            id="Top result",
        ),
        pytest.param([], {"text": None, "section": None}, id="No results"),
    ],
)
def test_infer_species(results: list[dict], expected: dict) -> None:
    """Test :func:`search.infer_species`.

    Args:
        results: Ranked search results.
        expected: Expected species chunk.

    """
    assert search.infer_species(results) == expected


@pytest.mark.parametrize(
    ("results", "expected"),
    [
        pytest.param(
            [
                {"section": "Conclusions", "text": "a", "score": 2.0},
                {"section": "Methods", "text": "b", "score": 1.0},
            ],
            {"text": "b", "section": "Methods", "score": 1.0},
            id="Excluded section skipped",
        ),
        pytest.param(
            [{"section": "Acknowledgements", "text": "a", "score": 2.0}],
            {"text": "a", "section": "Acknowledgements"},
            id="Only excluded sections, top result kept",
        ),
        pytest.param([], {"text": None, "section": None}, id="No results"),
    ],
)
def test_infer_strain(results: list[dict], expected: dict) -> None:
    """Test :func:`search.infer_strain`.

    Args:
        results: Ranked search results.
        expected: Expected strain chunk.

    """
    assert search.infer_strain(results) == expected


@pytest.mark.parametrize(
    ("rule", "vector", "expected"),
    [
        pytest.param(None, None, None, id="None"),
        pytest.param(None, "autopolyploid", "autopolyploid", id="One side"),
        pytest.param("allopolyploid", "amphidiploid", "amphidiploid", id="More specific wins"),
        pytest.param("unrecognised", "allopolyploid", "allopolyploid", id="Unknown term loses"),
    ],
)
def test_merge_mechanism(rule: str | None, vector: str | None, expected: str | None) -> None:
    """Test :func:`search._merge_mechanism`.

    Args:
        rule: Rule-based mechanism.
        vector: Vector-search mechanism.
        expected: Expected merged mechanism.

    """
    assert search._merge_mechanism(rule, vector) == expected


@pytest.mark.parametrize(
    ("level", "mechanism", "origin", "expected"),
    [
        pytest.param(
            4,
            "allopolyploid",
            {"status": "from_reference", "source": "goat"},
            (4, "allopolyploid", "tetraploid", "from_reference", "goat"),
            id="Status carried from origin",
        ),
        pytest.param(
            2,
            "allopolyploid",
            {},
            (2, None, "diploid", "found_in_paper", "paper"),
            id="Diploid, no mechanism",
        ),
        pytest.param(
            None,
            "autopolyploid",
            None,
            (None, "autopolyploid", None, "found_in_paper", "paper"),
            id="Mechanism only",
        ),
        pytest.param(None, None, {}, (None, None, None, "not_stated_in_paper", None), id="Nothing found"),
    ],
)
def test_build_ploidy_result(
    level: int | None, mechanism: str | None, origin: dict | None, expected: tuple
) -> None:
    """Test :func:`search._build_ploidy_result`.

    Args:
        level: Ploidy level.
        mechanism: Polyploidy mechanism.
        origin: Result the level came from.
        expected: Expected ``(level, mechanism, level_label, status, source)``.

    """
    result = search._build_ploidy_result(
        level, 0.7, "consensus", "evidence", origin, mechanism=mechanism, irregular=False, derivation="d"
    )
    fields = ("level", "mechanism", "level_label", "status", "source")
    assert tuple(result[field] for field in fields) == expected
    assert result["tuple"] == [result["level"], result["mechanism"], False, "d"]


@pytest.mark.parametrize(
    ("rule", "vector", "expected"),
    [
        pytest.param(
            {"mechanism": "allopolyploid", "confidence": 0.3, "evidence": "r"},
            {"irregular": True, "derivation": "x", "confidence": 0.2},
            (None, 0.4, "mechanism_only", "r", "allopolyploid", True, "x"),
            id="Mechanism only",
        ),
        pytest.param({}, {}, (None, 0.0, "both_unknown", None, None, False, None), id="Both unknown"),
        pytest.param(
            {"level": 4, "confidence": 0.6, "evidence": "r"},
            {"level": 4, "confidence": 0.8, "evidence": "v"},
            (4, 0.8, "consensus", "r", None, False, None),
            id="Consensus",
        ),
        pytest.param(
            {"level": 4, "confidence": 0.9},
            {"level": 4, "confidence": 0.95},
            (4, 0.95, "consensus", None, None, False, None),
            id="Consensus, capped",
        ),
        pytest.param(
            {},
            {"level": 6, "confidence": 0.5, "evidence": "v"},
            (6, 0.45, "vector_only", "v", None, False, None),
            id="Vector only",
        ),
        pytest.param(
            {"level": 6, "confidence": 0.5, "evidence": "r"},
            {},
            (6, 0.45, "rule_only", "r", None, False, None),
            id="Rule only",
        ),
        pytest.param(
            {"level": 2, "confidence": 0.5, "evidence": "r"},
            {"level": 4, "confidence": 0.5, "evidence": "v"},
            (2, 0.4, "rule_wins (conflict)", "r", None, False, None),
            id="Conflict, tie goes to rule",
        ),
        pytest.param(
            {"level": 2, "confidence": 0.5, "evidence": "r"},
            {"level": 4, "confidence": 0.6, "evidence": "v"},
            (4, 0.48, "vector_wins (conflict)", "v", None, False, None),
            id="Conflict, vector more confident",
        ),
    ],
)
def test_ensemble_ploidy_level(rule: dict, vector: dict, expected: tuple) -> None:
    """Test :func:`search._ensemble_ploidy_level`.

    Args:
        rule: Rule-based ploidy result.
        vector: Vector-search ploidy result.
        expected: Expected ``(level, confidence, method, evidence, mechanism, irregular, derivation)``.

    """
    result = search._ensemble_ploidy_level(rule, vector)
    fields = ("level", "confidence", "method", "evidence", "mechanism", "irregular", "derivation")
    assert tuple(result[field] for field in fields) == expected


@pytest.mark.parametrize(
    ("level", "rule", "vector", "expected"),
    [
        pytest.param(
            4,
            {
                "level": 4,
                "method": "chromosome_formula",
                "haploid_mentions": 1,
                "priority_derivation": "rule",
            },
            {"level": 4, "anchored": True, "haploid_mentions": 3, "priority_derivation": "vector"},
            (True, True, 3, "rule"),
            id="Formula-backed and anchored",
        ),
        pytest.param(
            4,
            {"level": 2, "method": "chromosome_formula", "anchored": True},
            {"priority_derivation": "vector"},
            (False, False, 0, "vector"),
            id="Support for another level",
        ),
        pytest.param(
            None,
            {"method": "chromosome_formula", "anchored": True},
            {},
            (False, False, 0, None),
            id="No level",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.search._ensemble_ploidy_level")
def test_ensemble_ploidy(
    mock_level: MagicMock, level: int | None, rule: dict, vector: dict, expected: tuple
) -> None:
    """Test :func:`search.ensemble_ploidy` sets the precision-guard flags for the chosen level.

    Args:
        mock_level: Patched `_ensemble_ploidy_level()`.
        level: Level chosen by the ensemble.
        rule: Rule-based ploidy result.
        vector: Vector-search ploidy result.
        expected: Expected ``(formula_backed, anchored, haploid_mentions, priority_derivation)``.

    """
    mock_level.return_value = {"level": level}
    result = search.ensemble_ploidy(rule, vector)
    fields = ("formula_backed", "anchored", "haploid_mentions", "priority_derivation")
    assert tuple(result[field] for field in fields) == expected


def test_ensemble_cultivars() -> None:
    """Test :func:`search.ensemble_cultivars`."""
    sections = {"Methods": "The cultivar Tifrunner was sequenced.", "Data availability": "Line 8-3"}
    cultivars = ["Tifrunner", "Line 8-3", "A cultivar name that is far too long to be real"]
    assert search.ensemble_cultivars(cultivars, {"section": "Methods"}, sections) == {
        "confirmed": ["Tifrunner"],
        "unconfirmed": ["Line 8-3"],
        "section": "Methods",
    }


@pytest.mark.parametrize(
    ("rule", "vector", "expected"),
    [
        pytest.param(
            {"value": "female"}, {"value": "female"}, {"value": "female", "method": "consensus"}, id="Agree"
        ),
        pytest.param({}, {"value": "male"}, {"value": "male", "method": "vector_only"}, id="Vector only"),
        pytest.param({"value": "male"}, {}, {"value": "male", "method": "rule_only"}, id="Rule only"),
        pytest.param(
            {"value": "male"}, {"value": "female"}, {"value": "unknown", "method": "conflict"}, id="Conflict"
        ),
    ],
)
def test_ensemble_sex(rule: dict, vector: dict, expected: dict) -> None:
    """Test :func:`search.ensemble_sex`.

    Args:
        rule: Rule-based sex call.
        vector: Vector-search sex call.
        expected: Expected combined call.

    """
    assert search.ensemble_sex(rule, vector) == expected


@pytest.mark.parametrize(
    ("rule_species", "species_text", "expected_species"),
    [
        pytest.param({"value": "Arachis hypogaea", "source": "ncbi"}, "x", None, id="Rule species kept"),
        pytest.param(None, "Arachis hypogaea L.", "Arachis hypogaea L.", id="Vector species fallback"),
    ],
)
@patch("ensembl.io.genomio.literature.search.infer_ploidy")
@patch("ensembl.io.genomio.literature.search.ncbi_chromosome_count")
@patch("ensembl.io.genomio.literature.search.search")
@patch("ensembl.io.genomio.literature.search.build_index")
def test_run_vector_search(
    mock_build: MagicMock,
    mock_search: MagicMock,
    mock_count: MagicMock,
    mock_infer_ploidy: MagicMock,
    rule_species: dict | None,
    species_text: str,
    expected_species: str | None,
) -> None:
    """Test :func:`search.run_vector_search` searches each field and ensembles with the rule result.

    Args:
        mock_build: Patched `build_index()`.
        mock_search: Patched `search()`.
        mock_count: Patched `ncbi_chromosome_count()`.
        mock_infer_ploidy: Patched `infer_ploidy()`.
        rule_species: Species from the rule-based stage.
        species_text: Text of the top species search result.
        expected_species: Expected vector-search species, or None when the rule species is kept.

    """
    mock_build.return_value = (np.zeros((1, 2)), MagicMock(), pd.DataFrame(), ["t"])
    ploidy_hits = [{"section": "Results", "text": f"p{i}", "score": float(i)} for i in range(6)]
    by_query = {
        "ploidy": ploidy_hits,
        "strain": [{"section": "Methods", "text": "Tifrunner", "score": 1.0}],
        "sex": [{"section": "Methods", "text": "a male plant", "score": 1.0}],
        "species": [{"section": "Results", "text": species_text, "score": 1.0}],
    }
    mock_search.side_effect = lambda query, *_, **__: by_query[
        next(k for k, q in search.QUERIES.items() if q == query)
    ]
    mock_count.return_value = 20
    mock_infer_ploidy.return_value = {"level": 4, "confidence": 0.8}
    rule = {
        "chromosome_number": {"value": 20, "source": "ncbi"},
        "target_names": ["arachis"],
        "ploidy": {"level": 4, "confidence": 0.6},
        "cultivars": ["Tifrunner"],
        "sex": {"value": "unknown"},
        "species": rule_species,
    }
    parsed = {"source": "PMC1", "chunks": [{"text": "t"}], "sections": {"Methods": "cultivar Tifrunner"}}

    final = search.run_vector_search(parsed, rule)

    mock_infer_ploidy.assert_called_once_with(ploidy_hits, chromosome_count=20, names=["arachis"])
    assert (final["ploidy"]["level"], final["ploidy"]["method"]) == (4, "consensus")
    assert final["ploidy"]["evidence_passages"] == [
        {"text": f"p{i}", "section": "Results", "score": float(i)} for i in range(5)
    ]
    assert final["cultivars"] == {"confirmed": ["Tifrunner"], "unconfirmed": [], "section": "Methods"}
    assert final["sex"] == {"value": "male", "method": "vector_only"}
    expected = rule_species or {"value": expected_species, "source": "vector_search"}
    assert final["species"] == expected
