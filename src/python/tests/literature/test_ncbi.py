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
"""Unit testing of `ensembl.io.genomio.literature.ncbi` module.

Typical usage example::
    $ pytest test_ncbi.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

from contextlib import nullcontext as does_not_raise
from typing import Any, ContextManager
from unittest.mock import MagicMock, call, patch

import pytest
import requests

from ensembl.io.genomio.literature import ncbi


def _http_error(status: int) -> MagicMock:
    """Build a fake response whose ``raise_for_status()`` raises an HTTP error.

    Args:
        status: HTTP status code to report.

    """
    response = MagicMock()
    response.raise_for_status.side_effect = requests.HTTPError(str(status))
    return response


@pytest.mark.parametrize(
    ("api_key", "side_effect", "expectation", "sleeps"),
    [
        pytest.param(None, [MagicMock()], does_not_raise(), [], id="First attempt"),
        pytest.param("KEY", [MagicMock()], does_not_raise(), [], id="API key added"),
        pytest.param(
            None, [requests.Timeout("slow"), MagicMock()], does_not_raise(), [call(1)], id="Retried once"
        ),
        pytest.param(
            None,
            [_http_error(503), MagicMock()],
            does_not_raise(),
            [call(1)],
            id="HTTP error status retried",
        ),
        pytest.param(
            None,
            [
                requests.ConnectionError("1"),
                requests.ConnectionError("2"),
                _http_error(500),
                _http_error(502),
            ],
            pytest.raises(requests.HTTPError, match="502"),
            [call(1), call(2), call(4)],
            id="All attempts fail, exponential backoff",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.ncbi.time.sleep")
@patch("ensembl.io.genomio.literature.ncbi.requests.get")
def test_entrez_get(
    mock_get: MagicMock,
    mock_sleep: MagicMock,
    api_key: str | None,
    side_effect: list,
    expectation: ContextManager,
    sleeps: list,
) -> None:
    """Test :func:`ncbi.entrez_get`.

    Args:
        mock_get: Patched `requests.get()`.
        mock_sleep: Patched `time.sleep()`, to check the backoff without waiting.
        api_key: Value of the module's NCBI API key.
        side_effect: Responses returned, or exceptions raised, by successive requests.
        expectation: Context manager for the expected outcome.
        sleeps: Expected backoff waits.

    """
    mock_get.side_effect = side_effect
    with patch.object(ncbi, "NCBI_API_KEY", api_key), expectation:
        assert ncbi.entrez_get("esearch.fcgi", {"db": "assembly"}, retries=4) is side_effect[-1]
    assert mock_get.call_args.kwargs["params"].get("api_key") == api_key
    assert mock_sleep.call_args_list == sleeps


@patch("ensembl.io.genomio.literature.ncbi.entrez_get")
def test_entrez_json(mock_entrez: MagicMock) -> None:
    """Test :func:`ncbi.entrez_json`.

    Args:
        mock_entrez: Patched `entrez_get()`.

    """
    mock_entrez.return_value.json.return_value = {"ok": True}
    assert ncbi.entrez_json("esearch.fcgi", {"db": "pmc"}) == {"ok": True}


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param([{"esearchresult": {"idlist": ["123", "456"]}}], "123", id="First UID"),
        pytest.param([{"esearchresult": {"idlist": []}}], None, id="No UID"),
        pytest.param(requests.ConnectionError("down"), None, id="Request error"),
    ],
)
@patch("ensembl.io.genomio.literature.ncbi.entrez_json")
def test_resolve_assembly_uid(mock_json: MagicMock, side_effect: Any, expected: str | None) -> None:
    """Test :func:`ncbi.resolve_assembly_uid`, including its per-accession cache.

    Args:
        mock_json: Patched `entrez_json()`.
        side_effect: Parsed esearch response, or the exception raised.
        expected: Expected assembly UID.

    """
    ncbi.resolve_assembly_uid.cache_clear()
    mock_json.side_effect = side_effect
    assert ncbi.resolve_assembly_uid("GCA_000001735.2") == expected
    assert ncbi.resolve_assembly_uid("GCA_000001735.2") == expected
    mock_json.assert_called_once()
    ncbi.resolve_assembly_uid.cache_clear()


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param(
            [
                {
                    "linksets": [
                        {
                            "linksetdbs": [
                                {"dbto": "nuccore", "links": ["9"]},
                                {"dbto": "pubmed", "links": ["1", "2"]},
                            ]
                        }
                    ]
                }
            ],
            ["1", "2"],
            id="Links to the requested database",
        ),
        pytest.param(
            [{"linksets": [{"linksetdbs": [{"dbto": "nuccore", "links": ["9"]}]}]}], [], id="No match"
        ),
        pytest.param([{"linksets": []}], [], id="Empty linksets"),
        pytest.param(requests.ConnectionError("down"), [], id="Request error"),
    ],
)
@patch("ensembl.io.genomio.literature.ncbi.entrez_json")
def test_entrez_elink(mock_json: MagicMock, side_effect: Any, expected: list) -> None:
    """Test :func:`ncbi._entrez_elink`.

    Args:
        mock_json: Patched `entrez_json()`.
        side_effect: Parsed elink response, or the exception raised.
        expected: Expected linked IDs.

    """
    mock_json.side_effect = side_effect
    assert ncbi._entrez_elink("assembly", "pubmed", "123") == expected


@pytest.mark.parametrize(
    ("uids", "side_effect", "expected"),
    [
        pytest.param([], None, {}, id="No UIDs, no request"),
        pytest.param(
            ["10", "20"],
            [
                MagicMock(
                    **{
                        "json.return_value": {
                            "linksets": [
                                {"ids": [10], "linksetdbs": [{"dbto": "pubmed", "links": ["1"]}]},
                                {"ids": [], "linksetdbs": [{"dbto": "pubmed", "links": ["x"]}]},
                                {"ids": [20], "linksetdbs": [{"dbto": "nuccore", "links": ["9"]}]},
                            ]
                        }
                    }
                )
            ],
            {"10": ["1"]},
            id="One linkset per UID",
        ),
        pytest.param(["10"], requests.ConnectionError("down"), {}, id="Request error"),
    ],
)
@patch("ensembl.io.genomio.literature.ncbi.entrez_get")
def test_entrez_elink_batch(mock_get: MagicMock, uids: list, side_effect: Any, expected: dict) -> None:
    """Test :func:`ncbi._entrez_elink_batch`.

    Args:
        mock_get: Patched `entrez_get()`.
        uids: UIDs to link from.
        side_effect: Elink response, or the exception raised.
        expected: Expected ``{uid: links}`` mapping.

    """
    mock_get.side_effect = side_effect
    assert ncbi._entrez_elink_batch("bioproject", "pubmed", uids) == expected


@pytest.mark.parametrize(
    ("uid", "pmids", "expected_titles"),
    [
        pytest.param(None, [], [], id="No assembly UID"),
        pytest.param("123", [], [], id="No linked PMIDs"),
        pytest.param("123", ["1", "2", "3"], ["Paper 1"], id="Missing paper skipped, max_results applied"),
    ],
)
@patch("ensembl.io.genomio.literature.ncbi.search_by_pmid")
@patch("ensembl.io.genomio.literature.ncbi._entrez_elink")
@patch("ensembl.io.genomio.literature.ncbi.resolve_assembly_uid")
def test_fetch_linked_pubmed_for_assembly(
    mock_uid: MagicMock,
    mock_elink: MagicMock,
    mock_pmid: MagicMock,
    uid: str | None,
    pmids: list,
    expected_titles: list,
) -> None:
    """Test :func:`ncbi.fetch_linked_pubmed_for_assembly`.

    Args:
        mock_uid: Patched `resolve_assembly_uid()`.
        mock_elink: Patched `_entrez_elink()`.
        mock_pmid: Patched `search_by_pmid()`; PMID "2" has no Europe PMC record.
        uid: Assembly UID returned by the lookup.
        pmids: PMIDs linked to the assembly.
        expected_titles: Titles of the expected papers.

    """
    mock_uid.return_value = uid
    mock_elink.return_value = pmids
    mock_pmid.side_effect = lambda pmid: None if pmid == "2" else {"title": f"Paper {pmid}"}
    papers = ncbi.fetch_linked_pubmed_for_assembly("GCA_000001735.2", max_results=2)
    assert [paper["title"] for paper in papers] == expected_titles


@pytest.mark.parametrize(
    ("uid", "bioprojects", "links", "expected_pmids"),
    [
        pytest.param(None, [], {}, [], id="No assembly UID"),
        pytest.param("123", [], {}, [], id="No BioProject"),
        pytest.param("123", ["7"], {}, [], id="No reference PMID"),
        pytest.param(
            "123",
            ["7", "8", "9", "10"],
            {"7": ["1"], "8": ["1", "2", "3"], "10": ["4"]},
            ["1", "2"],
            id="PMIDs de-duplicated, max_results applied",
        ),
        pytest.param("123", ["7", "8", "9", "10"], {"10": ["4"]}, [], id="Only the first three BioProjects"),
    ],
)
@patch("ensembl.io.genomio.literature.ncbi.search_by_pmid")
@patch("ensembl.io.genomio.literature.ncbi._entrez_elink_batch")
@patch("ensembl.io.genomio.literature.ncbi._entrez_elink")
@patch("ensembl.io.genomio.literature.ncbi.resolve_assembly_uid")
def test_fetch_bioproject_reference_papers(
    mock_uid: MagicMock,
    mock_elink: MagicMock,
    mock_batch: MagicMock,
    mock_pmid: MagicMock,
    uid: str | None,
    bioprojects: list,
    links: dict,
    expected_pmids: list,
) -> None:
    """Test :func:`ncbi.fetch_bioproject_reference_papers`.

    Args:
        mock_uid: Patched `resolve_assembly_uid()`.
        mock_elink: Patched `_entrez_elink()`.
        mock_batch: Patched `_entrez_elink_batch()`.
        mock_pmid: Patched `search_by_pmid()`.
        uid: Assembly UID returned by the lookup.
        bioprojects: BioProject UIDs linked to the assembly.
        links: Reference PMIDs per BioProject.
        expected_pmids: PMIDs of the expected papers, in order.

    """
    mock_uid.return_value = uid
    mock_elink.return_value = bioprojects
    mock_batch.side_effect = lambda _from, _to, uids: {uid: links[uid] for uid in uids if uid in links}
    mock_pmid.side_effect = lambda pmid: {"pmid": pmid}
    papers = ncbi.fetch_bioproject_reference_papers("GCA_000001735.2", max_results=2)
    assert papers == [
        {"pmid": pmid, "retrieval_source": "bioproject", "is_reference_paper": True}
        for pmid in expected_pmids
    ]


# Shape of the "articleids" list in a real PMC esummary record (PMC6175935).
ARTICLE_IDS = [
    {"idtype": "pmid", "value": "30297696"},
    {"idtype": "pmcid", "value": "PMC6175935"},
    {"idtype": "doi", "value": "10.1038/s41598-018-33344-8"},
]


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param(requests.ConnectionError("down"), [], id="Search error"),
        pytest.param([{"esearchresult": {"idlist": []}}], [], id="No PMC IDs"),
        pytest.param(
            [{"esearchresult": {"idlist": ["1", "2"]}}, requests.ConnectionError("down")],
            [],
            id="Summary error",
        ),
        pytest.param(
            [
                {"esearchresult": {"idlist": ["1", "2"]}},
                {"result": {"1": {"title": "Peanut genome", "articleids": ARTICLE_IDS}, "2": {}}},
            ],
            [
                {
                    "title": "Peanut genome",
                    "pmcid": "PMC1",
                    "pmid": "30297696",
                    "doi": "10.1038/s41598-018-33344-8",
                    "abstractText": "",
                    "source": "ncbi_entrez",
                }
            ],
            id="Summaries converted, IDs read from articleids, empty skipped",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.ncbi.entrez_json")
def test_search_ncbi_entrez(mock_json: MagicMock, side_effect: Any, expected: list) -> None:
    """Test :func:`ncbi.search_ncbi_entrez`.

    Args:
        mock_json: Patched `entrez_json()`.
        side_effect: Parsed esearch and esummary responses, or the exception raised.
        expected: Expected paper records.

    """
    mock_json.side_effect = side_effect
    assert ncbi.search_ncbi_entrez("Arachis hypogaea", max_results=2) == expected
    term = mock_json.call_args_list[0].args[1]["term"]
    assert term == '"Arachis hypogaea" AND "genome assembly" AND "open access"[filter]'
