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
"""Unit testing of `ensembl.io.genomio.literature.fetch` module.

Typical usage example::
    $ pytest test_fetch.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

from typing import Any
from unittest.mock import MagicMock, call, patch

import pytest
import requests

from ensembl.io.genomio.literature import fetch

PEANUT = {"scientific_name": "Arachis hypogaea", "common_name": "peanut"}
ACCESSION = "GCA_003086295.3"


@pytest.mark.parametrize(
    ("paper", "identifiers", "expected"),
    [
        pytest.param({}, PEANUT, 0.0, id="Nothing"),
        pytest.param(
            {"title": "Genome assembly of Arachis hypogaea"},
            PEANUT,
            6.0,
            id="Name and genome keyword in title",
        ),
        pytest.param(
            {"title": "Assembling the genome", "abstractText": "Arachis hypogaea"},
            PEANUT,
            4.0,
            id="Name in abstract, genome and assembly words in title",
        ),
        pytest.param({"title": "Peanut genome"}, PEANUT, 2.5, id="Common name and genome word in title"),
        pytest.param(
            {"title": "Pea genome"},
            {"scientific_name": "Pisum sativum", "common_name": "pea"},
            1.5,
            id="Common name too short",
        ),
        pytest.param(
            {"title": "A study", "abstractText": "whole genome of a tetraploid"},
            PEANUT,
            2.0,
            id="Genome keyword and ploidy hint in abstract",
        ),
        pytest.param({"pmCid": "PMC1"}, PEANUT, 2.0, id="Open access"),
        pytest.param(
            {"retrieval_source": "elink", "abstractText": "Arachis hypogaea"},
            PEANUT,
            6.0,
            id="Accession-linked and names the species",
        ),
        pytest.param({"is_reference_paper": True}, PEANUT, 2.5, id="BioProject-linked"),
        pytest.param({"_fulltext_name_confirmed": True}, PEANUT, 4.0, id="Full-text confirmed"),
        pytest.param({"title": "Chloroplast genome"}, PEANUT, -8.5, id="Organelle in title"),
        pytest.param({"abstractText": "plastome"}, PEANUT, -2.0, id="Organelle in abstract"),
        pytest.param({"title": "Bacillus isolated from peanut"}, PEANUT, -9.0, id="Host context in title"),
        pytest.param({"abstractText": "rhizosphere"}, PEANUT, -2.0, id="Host context in abstract"),
    ],
)
def test_score_paper_candidate(paper: dict, identifiers: dict, expected: float) -> None:
    """Test :func:`fetch.score_paper_candidate`.

    Args:
        paper: Candidate paper record.
        identifiers: Assembly identifiers.
        expected: Expected relevance score.

    """
    assert fetch.score_paper_candidate(paper, identifiers) == expected


@patch("ensembl.io.genomio.literature.fetch.score_paper_candidate")
def test_dedupe_and_rank(mock_score: MagicMock) -> None:
    """Test :func:`fetch._dedupe_and_rank`.

    Args:
        mock_score: Patched `score_paper_candidate()`, returning each paper's ``s``.

    """
    mock_score.side_effect = lambda paper, _: paper["s"]
    candidates = [
        {"id": "a", "pmid": "1", "s": 1.0},
        {"id": "b", "pmid": " 1 ", "s": 9.0},
        {"id": "c", "pmCid": "PMC2", "s": 3.333},
        {"id": "d", "pmcid": "PMC2", "s": 9.0},
        {"id": "e", "doi": "10.1/x", "s": 2.0},
        {"id": "f", "title": "A Long Title", "s": 5.0},
        {"id": "g", "title": "a long title", "s": 9.0},
        {"id": "h", "s": 9.0},
    ]
    ranked = fetch._dedupe_and_rank(candidates, PEANUT)
    assert [(paper["id"], paper["relevance_score"]) for paper in ranked] == [
        ("f", 5.0),
        ("c", 3.33),
        ("e", 2.0),
        ("a", 1.0),
    ]


@pytest.mark.parametrize(
    ("papers", "expected"),
    [
        pytest.param([{"title": "Genome"}], False, id="No PMCID"),
        pytest.param([{"pmcid": "PMC1", "title": "Mitochondria genome"}], False, id="Organelle only"),
        pytest.param(
            [{"pmcid": "PMC1", "abstractText": "chloroplast"}, {"pmCid": "PMC2", "title": "Genome"}],
            True,
            id="One usable",
        ),
    ],
)
def test_has_usable_fulltext(papers: list[dict], expected: bool) -> None:
    """Test :func:`fetch.has_usable_fulltext`.

    Args:
        papers: Candidate papers.
        expected: Whether a usable open-access nuclear-genome paper is present.

    """
    assert fetch.has_usable_fulltext(papers) == expected


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param(
            [
                MagicMock(
                    **{
                        "json.return_value": {
                            "data": [
                                {
                                    "title": "Arachis genome",
                                    "abstract": None,
                                    "externalIds": {"PubMedCentral": "1", "PubMed": 2, "DOI": "10.1/x"},
                                },
                                {"title": "Unrelated", "abstract": "Oryza", "externalIds": {}},
                                {"title": "Second", "abstract": "arachis", "externalIds": {}},
                                {"title": "Third", "abstract": "Arachis", "externalIds": {}},
                            ]
                        }
                    }
                )
            ],
            [
                {
                    "title": "Arachis genome",
                    "pmcid": "PMC1",
                    "pmid": "2",
                    "doi": "10.1/x",
                    "abstractText": "",
                    "source": "semantic_scholar",
                },
                {
                    "title": "Second",
                    "pmcid": "",
                    "pmid": "",
                    "doi": "",
                    "abstractText": "arachis",
                    "source": "semantic_scholar",
                },
            ],
            id="Genus filter, IDs mapped, capped",
        ),
        pytest.param(requests.ConnectionError("down"), [], id="Request error"),
    ],
)
@patch("ensembl.io.genomio.literature.fetch.requests.get")
def test_search_semantic_scholar(mock_get: MagicMock, side_effect: Any, expected: list) -> None:
    """Test :func:`fetch.search_semantic_scholar`.

    Args:
        mock_get: Patched `requests.get()`.
        side_effect: Response returned, or exception raised, by the request.
        expected: Expected paper records.

    """
    mock_get.side_effect = side_effect
    assert fetch.search_semantic_scholar("Arachis hypogaea", max_results=2) == expected


@pytest.mark.xfail(strict=True, raises=IndexError, reason="fetch.py:252 indexes ''.split()[0]")
@patch("ensembl.io.genomio.literature.fetch.requests.get")
def test_search_semantic_scholar_without_name(mock_get: MagicMock) -> None:
    """Test :func:`fetch.search_semantic_scholar` returns [] without a scientific name, as documented.

    Reachable from `fetch_papers_for_assembly` when an assembly has a taxon ID but no scientific name.

    Args:
        mock_get: Patched `requests.get()`.

    """
    mock_get.return_value.json.return_value = {"data": [{"title": "Any paper", "externalIds": {}}]}
    assert not fetch.search_semantic_scholar("")


@patch("ensembl.io.genomio.literature.fetch.score_paper_candidate")
def test_strongest_score(mock_score: MagicMock) -> None:
    """Test :func:`fetch._strongest_score`.

    Args:
        mock_score: Patched `score_paper_candidate()`, returning each paper's ``s``.

    """
    mock_score.side_effect = lambda paper, _: paper["s"]
    assert fetch._strongest_score([], PEANUT) == 0.0
    assert fetch._strongest_score([{"s": -2.0}, {"s": 4.5}, {"s": 1.0}], PEANUT) == 4.5


@pytest.mark.parametrize(
    ("paper", "names", "expected"),
    [
        pytest.param({}, [], True, id="No names to check"),
        pytest.param({"title": "The Arachis genome"}, ["arachis"], True, id="Name in title"),
        pytest.param({"text_data": {"text": "...peanuts are..."}}, ["peanut"], True, id="Word start in text"),
        pytest.param({"title": "Price of rice"}, ["ice"], False, id="Not at a word start"),
        pytest.param(
            {"text_data": {"text": "x" * 3000 + " arachis"}}, ["arachis"], False, id="Beyond the head"
        ),
    ],
)
def test_names_target_up_front(paper: dict, names: list[str], expected: bool) -> None:
    """Test :func:`fetch._names_target_up_front`.

    Args:
        paper: Paper record.
        names: Target organism names.
        expected: Whether the paper names the target up front.

    """
    assert fetch._names_target_up_front(paper, names) == expected


MANY_OTHERS = [f"gca_00000000{i}.1" for i in range(5)]


@pytest.mark.parametrize(
    ("paper", "expected", "others"),
    [
        pytest.param({"retrieval_source": "linked_pmid"}, True, None, id="Linked PMID"),
        pytest.param({"retrieval_source": "elink"}, True, None, id="Entrez elink"),
        pytest.param(
            {"retrieval_source": "bioproject", "title": "Arachis genome"}, True, None, id="BioProject"
        ),
        pytest.param(
            {"is_reference_paper": True, "title": "RefSeq"}, False, None, id="BioProject, other topic"
        ),
        pytest.param(
            {"title": "A", "text_data": {"text": "Data: gca_003086295 and GCF_003086295.1, GCA_1.2"}},
            True,
            ["gca_1.2", "gcf_003086295.1"],
            id="Accession in text, other accessions recorded",
        ),
        pytest.param({"abstractText": "GCA_9.1"}, False, ["gca_9.1"], id="Accession absent"),
        pytest.param(
            {"title": "Many plants", "text_data": {"text": " ".join(["GCA_003086295"] + MANY_OTHERS)}},
            False,
            MANY_OTHERS,
            id="Among many others, organism not named",
        ),
        pytest.param(
            {"title": "Arachis and others", "text_data": {"text": " ".join(["GCA_003086295"] + MANY_OTHERS)}},
            True,
            MANY_OTHERS,
            id="Among many others, organism named",
        ),
    ],
)
def test_paper_is_accession_linked(paper: dict, expected: bool, others: list[str] | None) -> None:
    """Test :func:`fetch._paper_is_accession_linked`.

    Args:
        paper: Prepared paper record.
        expected: Whether the paper is accepted.
        others: Expected ``other_accessions_in_text``, if any are recorded.

    """
    assert fetch._paper_is_accession_linked(paper, ACCESSION, ["arachis", "peanut"]) == expected
    assert paper.get("other_accessions_in_text") == others


@pytest.mark.parametrize(
    ("linked_pmids", "sources", "expected", "called"),
    [
        pytest.param(
            ["1", "2"],
            {"pmid": [{"id": "p1", "s": 8.0}, None], "bioproject": [{"id": "bp", "s": 0.0}]},
            [("p1", "linked_pmid"), ("bp", None)],
            ["pmid", "bioproject"],
            id="Strong linked paper stops escalation",
        ),
        pytest.param(
            [],
            {"elink": [{"id": "e1", "s": 6.0}], "fulltext": [{"id": "f1", "s": 0.0}]},
            [("e1", "elink"), ("f1", "europepmc_fulltext")],
            ["elink", "bioproject", "fulltext"],
            id="Medium paper adds full-text search only",
        ),
        pytest.param(
            [],
            {
                "name": [{"id": "n1", "s": 1.0, "retrieval_source": "kept"}],
                "entrez": [{"id": "x1", "s": 0.0}],
            },
            [("n1", "kept"), ("x1", "entrez")],
            ["elink", "bioproject", "name", "fulltext", "entrez"],
            id="Weak paper escalates to Entrez",
        ),
        pytest.param(
            [],
            {"semantic": [{"id": "s1", "s": 0.0}]},
            [("s1", "semantic_scholar")],
            ["elink", "bioproject", "name", "fulltext", "entrez", "semantic"],
            id="Nothing found escalates to Semantic Scholar",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.fetch.score_paper_candidate")
@patch("ensembl.io.genomio.literature.fetch.search_semantic_scholar")
@patch("ensembl.io.genomio.literature.fetch.search_ncbi_entrez")
@patch("ensembl.io.genomio.literature.fetch.search_europe_pmc_fulltext")
@patch("ensembl.io.genomio.literature.fetch.search_europe_pmc_by_name")
@patch("ensembl.io.genomio.literature.fetch.fetch_bioproject_reference_papers")
@patch("ensembl.io.genomio.literature.fetch.fetch_linked_pubmed_for_assembly")
@patch("ensembl.io.genomio.literature.fetch.search_by_pmid")
def test_gather_paper_candidates(
    mock_pmid: MagicMock,
    mock_elink: MagicMock,
    mock_bioproject: MagicMock,
    mock_name: MagicMock,
    mock_fulltext: MagicMock,
    mock_entrez: MagicMock,
    mock_semantic: MagicMock,
    mock_score: MagicMock,
    linked_pmids: list[str],
    sources: dict[str, list],
    expected: list[tuple],
    called: list[str],
) -> None:
    """Test :func:`fetch._gather_paper_candidates` escalates through the sources in priority order.

    Args:
        mock_pmid: Patched `search_by_pmid()`.
        mock_elink: Patched `fetch_linked_pubmed_for_assembly()`.
        mock_bioproject: Patched `fetch_bioproject_reference_papers()`.
        mock_name: Patched `search_europe_pmc_by_name()`.
        mock_fulltext: Patched `search_europe_pmc_fulltext()`.
        mock_entrez: Patched `search_ncbi_entrez()`.
        mock_semantic: Patched `search_semantic_scholar()`.
        mock_score: Patched `score_paper_candidate()`, returning each paper's ``s``.
        linked_pmids: PMIDs linked from the assembly record.
        sources: Papers returned by each source (the ``pmid`` source gives one result per PMID).
        expected: Expected ``(id, retrieval_source)`` of each candidate, in order.
        called: Sources expected to be queried.

    """
    mocks = {
        "pmid": mock_pmid,
        "elink": mock_elink,
        "bioproject": mock_bioproject,
        "name": mock_name,
        "fulltext": mock_fulltext,
        "entrez": mock_entrez,
        "semantic": mock_semantic,
    }
    for name, mock in mocks.items():
        mock.side_effect = sources["pmid"] if name == "pmid" and "pmid" in sources else None
        mock.return_value = sources.get(name, [])
    mock_score.side_effect = lambda paper, _: paper["s"]
    identifiers = {"scientific_name": "Arachis hypogaea", "common_name": "peanut", "taxon_id": "3818"}
    candidates = fetch._gather_paper_candidates(ACCESSION, {"linked_pmids": linked_pmids}, identifiers, 3)
    assert [(paper["id"], paper.get("retrieval_source")) for paper in candidates] == expected
    assert [name for name, mock in mocks.items() if mock.called] == called


@patch("ensembl.io.genomio.literature.fetch.get_fulltext_by_pmcid")
def test_confirm_fulltext_candidates(mock_fulltext: MagicMock) -> None:
    """Test :func:`fetch._confirm_fulltext_candidates`.

    Args:
        mock_fulltext: Patched `get_fulltext_by_pmcid()`; only PMC1 and PMC6 mention the species.

    """
    mock_fulltext.side_effect = lambda pmcid: "... ARACHIS HYPOGAEA ..." if pmcid in ("PMC1", "PMC6") else "x"
    candidates = [
        {"id": "plain"},
        {"id": "1", "pmcid": "PMC1", "_needs_fulltext_confirm": True},
        {"id": "no-pmcid", "_needs_fulltext_confirm": True},
    ] + [{"id": str(i), "pmCid": f"PMC{i}", "_needs_fulltext_confirm": True} for i in range(2, 7)]
    confirmed = fetch._confirm_fulltext_candidates(candidates, "Arachis hypogaea")
    assert confirmed == [{"id": "plain"}, {"id": "1", "pmcid": "PMC1", "_fulltext_name_confirmed": True}]
    assert mock_fulltext.call_args_list == [call(f"PMC{i}") for i in range(1, 6)]


@patch("ensembl.io.genomio.literature.fetch.get_best_available_text")
@patch("ensembl.io.genomio.literature.fetch.validate_paper_ids")
def test_prepare_papers(mock_validate: MagicMock, mock_text: MagicMock) -> None:
    """Test :func:`fetch._prepare_papers`.

    Args:
        mock_validate: Patched `validate_paper_ids()`.
        mock_text: Patched `get_best_available_text()`.

    """
    mock_validate.return_value = {"best_available": "fulltext"}
    mock_text.return_value = {"source": "fulltext", "text": "<article/>"}
    paper = {"title": "T", "pmcid": "PMC1", "retrieval_source": "elink", "relevance_score": 7.5, "extra": 1}
    assert fetch._prepare_papers([paper]) == [
        {
            "title": "T",
            "pmcid": "PMC1",
            "pmid": None,
            "doi": None,
            "retrieval_source": "elink",
            "relevance_score": 7.5,
            "is_reference_paper": False,
            "validation": {"best_available": "fulltext"},
            "text_data": {"source": "fulltext", "text": "<article/>"},
        }
    ]


@pytest.mark.parametrize(
    ("papers", "expected"),
    [
        pytest.param([{"title": "Name hit"}], None, id="None"),
        pytest.param(
            [
                {"title": "Name hit"},
                {"title": "Ref", "pmid": "1", "is_reference_paper": True, "relevance_score": 4},
            ],
            {
                "title": "Ref",
                "pmcid": None,
                "pmid": "1",
                "doi": None,
                "retrieval_source": "bioproject",
                "relevance_score": 4,
            },
            id="Reference paper",
        ),
    ],
)
def test_extract_reference_paper(papers: list[dict], expected: dict | None) -> None:
    """Test :func:`fetch._extract_reference_paper`.

    Args:
        papers: Ranked papers.
        expected: Expected reference paper record.

    """
    assert fetch._extract_reference_paper(papers) == expected


PREPARED = [{"id": "keep", "ok": True}, {"id": "drop", "ok": False}]


@pytest.mark.parametrize(
    ("assembly", "expected"),
    [
        pytest.param(
            {"assembly_accession": ACCESSION, "scientific_name": "Arachis hypogaea", "taxon_id": "3818"},
            {"papers": [{"id": "keep", "ok": True}], "reference_paper": {"title": "Ref"}},
            id="Only accession-linked papers kept",
        ),
        pytest.param({"assembly_accession": ACCESSION}, {"papers": []}, id="No name or taxon, no search"),
    ],
)
@patch(
    "ensembl.io.genomio.literature.fetch._extract_reference_paper", MagicMock(return_value={"title": "Ref"})
)
@patch("ensembl.io.genomio.literature.fetch._paper_is_accession_linked", lambda paper, *_: paper["ok"])
@patch("ensembl.io.genomio.literature.fetch._prepare_papers", MagicMock(return_value=PREPARED))
@patch("ensembl.io.genomio.literature.fetch._dedupe_and_rank", MagicMock())
@patch("ensembl.io.genomio.literature.fetch._confirm_fulltext_candidates", MagicMock())
@patch("ensembl.io.genomio.literature.fetch._gather_paper_candidates", MagicMock())
@patch("ensembl.io.genomio.literature.fetch.fetch_assembly_metadata", MagicMock())
@patch("ensembl.io.genomio.literature.fetch.enrich_assembly_metadata")
def test_fetch_papers_for_assembly(mock_enrich: MagicMock, assembly: dict, expected: dict) -> None:
    """Test :func:`fetch.fetch_papers_for_assembly` keeps only accession-linked papers.

    Args:
        mock_enrich: Patched `enrich_assembly_metadata()`.
        assembly: Enriched assembly metadata.
        expected: Expected papers (and reference paper).

    """
    mock_enrich.return_value = assembly
    result = fetch.fetch_papers_for_assembly(ACCESSION)
    assert {"assembly": assembly, **expected}.items() <= result.items()


@pytest.mark.xfail(strict=True, reason="fetch.py:539 early return omits the 'reference_paper' key")
@patch("ensembl.io.genomio.literature.fetch.enrich_assembly_metadata", lambda assembly: assembly)
@patch("ensembl.io.genomio.literature.fetch.fetch_assembly_metadata", lambda accession: {})
def test_fetch_papers_for_assembly_shape_without_names() -> None:
    """Test :func:`fetch.fetch_papers_for_assembly` always returns the documented keys."""
    assert set(fetch.fetch_papers_for_assembly(ACCESSION)) == {"assembly", "papers", "reference_paper"}
