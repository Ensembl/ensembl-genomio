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
"""Unit testing of `ensembl.io.genomio.literature.europepmc` module.

Typical usage example::
    $ pytest test_europepmc.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

import io
from typing import Any
from unittest.mock import MagicMock, call, patch
import zipfile

import pytest
import requests

from ensembl.io.genomio.literature import europepmc


def _response(payload: Any = None, status: int = 200, text: str = "", content: bytes = b"") -> MagicMock:
    """Build a fake `requests.Response`.

    Args:
        payload: Value returned by ``json()``.
        status: HTTP status code.
        text: Response body as text.
        content: Response body as bytes.

    """
    response = MagicMock(status_code=status, text=text, content=content)
    response.json.return_value = payload
    if status >= 500:
        response.raise_for_status.side_effect = requests.HTTPError(str(status))
    return response


def _zip(members: dict[str, bytes]) -> bytes:
    """Build an in-memory zip archive.

    Args:
        members: File name to content.

    """
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, data in members.items():
            archive.writestr(name, data)
    return buffer.getvalue()


@pytest.mark.parametrize(
    ("side_effect", "expected", "sleeps"),
    [
        pytest.param([_response({"resultList": {"result": [{"id": "1"}]}})], [{"id": "1"}], [], id="Results"),
        pytest.param([_response({})], [], [], id="No result list"),
        pytest.param(
            [
                requests.Timeout("slow"),
                _response(status=503),
                _response({"resultList": {"result": [{"id": "1"}]}}),
            ],
            [{"id": "1"}],
            [call(1), call(2)],
            id="Timeout and HTTP error retried",
        ),
        pytest.param([requests.ConnectionError("down")] * 4, [], [call(1), call(2), call(4)], id="All fail"),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.time.sleep")
@patch("ensembl.io.genomio.literature.europepmc.requests.get")
def test_search_europe_pmc(
    mock_get: MagicMock, mock_sleep: MagicMock, side_effect: list, expected: list, sleeps: list
) -> None:
    """Test :func:`europepmc.search_europe_pmc`.

    Args:
        mock_get: Patched `requests.get()`.
        mock_sleep: Patched `time.sleep()`, to check the backoff without waiting.
        side_effect: Responses returned, or exceptions raised, by successive requests.
        expected: Expected paper records.
        sleeps: Expected backoff waits.

    """
    mock_get.side_effect = side_effect
    assert europepmc.search_europe_pmc("wheat", retries=4) == expected
    assert mock_sleep.call_args_list == sleeps


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param(
            [_response({"resultList": {"result": [{"pmid": "1"}, {"pmid": "2"}]}})], {"pmid": "1"}, id="Found"
        ),
        pytest.param([_response({"resultList": {"result": []}})], None, id="Not found"),
        pytest.param(requests.ConnectionError("down"), None, id="Request error"),
        pytest.param([MagicMock(**{"json.side_effect": ValueError("bad json")})], None, id="Invalid JSON"),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.requests.get")
def test_search_by_pmid(mock_get: MagicMock, side_effect: Any, expected: dict | None) -> None:
    """Test :func:`europepmc.search_by_pmid`.

    Args:
        mock_get: Patched `requests.get()`.
        side_effect: Response returned, or exception raised, by the request.
        expected: Expected paper record.

    """
    mock_get.side_effect = side_effect
    assert europepmc.search_by_pmid("1") == expected
    assert mock_get.call_args.kwargs["params"]["query"] == "EXT_ID:1 AND SRC:MED"


@pytest.mark.parametrize(
    ("names", "results", "expected", "queries"),
    [
        pytest.param(
            ("3818", "Arachis hypogaea", "peanut"),
            [
                [
                    {"title": "Arachis hypogaea genome"},
                    {"title": "Peanut allergy"},
                    {"abstractText": "arachis"},
                    {"title": "Arachis duranensis"},
                ]
            ],
            [{"title": "Arachis hypogaea genome"}, {"abstractText": "arachis"}],
            ['("Arachis hypogaea" OR "peanut") AND "genome assembly"'],
            id="First query, genus filter, capped",
        ),
        pytest.param(
            ("", "Drosophila melanogaster", "fly"),
            [[{"title": "Fly behaviour"}], [{"title": "Drosophila genome"}]],
            [{"title": "Drosophila genome"}],
            [
                '"Drosophila melanogaster" AND "genome assembly"',
                '"Drosophila melanogaster" AND ("genome sequence" OR ploidy OR chromosome)',
            ],
            id="Short common name ignored, irrelevant hits escalate",
        ),
        pytest.param(
            ("4530", "Oryza sativa var. japonica", ""),
            [[]] * 6,
            [],
            [
                '"Oryza sativa var. japonica" AND "genome assembly"',
                '"Oryza sativa var. japonica" AND ("genome sequence" OR ploidy OR chromosome)',
                '"Oryza sativa var. japonica" AND "genome sequence"',
                '"Oryza sativa" AND "genome assembly"',
                '"Oryza sativa" AND ("genome sequence" OR ploidy OR chromosome)',
                "TAXONOMY:4530 AND genome assembly",
            ],
            id="Trinomial falls back to binomial, then taxon ID",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.search_europe_pmc")
def test_search_europe_pmc_by_name(
    mock_search: MagicMock, names: tuple[str, str, str], results: list, expected: list, queries: list[str]
) -> None:
    """Test :func:`europepmc.search_europe_pmc_by_name`.

    Args:
        mock_search: Patched `search_europe_pmc()`.
        names: Taxon ID, scientific name and common name.
        results: Results returned by each successive query.
        expected: Expected paper records.
        queries: Queries expected to be sent, in order.

    """
    mock_search.side_effect = results
    taxon_id, scientific_name, common_name = names
    assert (
        europepmc.search_europe_pmc_by_name(taxon_id, scientific_name, common_name, max_results=2) == expected
    )
    assert mock_search.call_args_list == [call(query, max_results=6) for query in queries]


@pytest.mark.parametrize(
    ("name", "results", "expected_ids"),
    [
        pytest.param("", [], [], id="No name"),
        pytest.param(
            "Agapanthus africanus",
            [
                [
                    {"id": "a", "pmcid": "PMC1"},
                    {"id": "b"},
                    {"id": "c", "pmCid": "PMC3"},
                    {"id": "d", "pmcid": "4"},
                ]
            ],
            ["a", "c"],
            id="Open access only, capped",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.search_europe_pmc")
def test_search_europe_pmc_fulltext(
    mock_search: MagicMock, name: str, results: list, expected_ids: list[str]
) -> None:
    """Test :func:`europepmc.search_europe_pmc_fulltext`.

    Args:
        mock_search: Patched `search_europe_pmc()`.
        name: Scientific name to search for.
        results: Results returned by the search.
        expected_ids: IDs of the expected candidates.

    """
    mock_search.side_effect = results
    candidates = europepmc.search_europe_pmc_fulltext(name, max_results=1)
    assert [paper["id"] for paper in candidates] == expected_ids
    assert all(paper["_needs_fulltext_confirm"] for paper in candidates)


@pytest.mark.parametrize(
    ("paper", "head", "pmid_record", "expected"),
    [
        pytest.param(
            {"pmcid": "PMC1", "pmid": "1", "doi": "10.1/x"},
            _response(status=200),
            {"pmid": "1"},
            (True, True, True, "fulltext"),
            id="Full text",
        ),
        pytest.param(
            {"pmcid": "PMC1", "pmid": "1"},
            _response(status=404),
            {"pmid": "1"},
            (False, True, False, "abstract"),
            id="Abstract only",
        ),
        pytest.param(
            {"pmcid": "PMC1", "pmid": "1"},
            requests.ConnectionError("down"),
            None,
            (False, False, False, "none"),
            id="Nothing retrievable",
        ),
        pytest.param({}, None, None, (False, False, False, "none"), id="No IDs"),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.search_by_pmid")
@patch("ensembl.io.genomio.literature.europepmc.requests.head")
def test_validate_paper_ids(
    mock_head: MagicMock,
    mock_pmid: MagicMock,
    paper: dict,
    head: Any,
    pmid_record: dict | None,
    expected: tuple,
) -> None:
    """Test :func:`europepmc.validate_paper_ids`.

    Args:
        mock_head: Patched `requests.head()` for the full-text check.
        mock_pmid: Patched `search_by_pmid()`.
        paper: Paper record to validate.
        head: Response returned, or exception raised, by the full-text check.
        pmid_record: Record returned by the PMID lookup.
        expected: Expected ``(pmcid_valid, pmid_valid, doi_present, best_available)``.

    """
    mock_head.side_effect = [head]
    mock_pmid.return_value = pmid_record
    result = europepmc.validate_paper_ids(paper)
    fields = ("pmcid_valid", "pmid_valid", "doi_present", "best_available")
    assert tuple(result[field] for field in fields) == expected


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param([_response(status=200, text="<article/>")], "<article/>", id="Found"),
        pytest.param([_response(status=404, text="Not found")], None, id="Not found"),
        pytest.param(requests.ConnectionError("down"), None, id="Request error"),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.requests.get")
def test_get_fulltext_by_pmcid(mock_get: MagicMock, side_effect: Any, expected: str | None) -> None:
    """Test :func:`europepmc.get_fulltext_by_pmcid`.

    Args:
        mock_get: Patched `requests.get()`.
        side_effect: Response returned, or exception raised, by the request.
        expected: Expected full-text XML.

    """
    mock_get.side_effect = side_effect
    assert europepmc.get_fulltext_by_pmcid("PMC1") == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        pytest.param(_zip({"word/document.xml": b"<w:t>2n</w:t>"}), "<w:t>2n</w:t>", id="Member read"),
        pytest.param(_zip({"other.xml": b"x"}), "", id="Member missing"),
        pytest.param(b"not a zip", "", id="Not a zip"),
    ],
)
def test_read_zip(raw: bytes, expected: str) -> None:
    """Test :func:`europepmc._read_zip`.

    Args:
        raw: Archive bytes.
        expected: Expected member text.

    """
    assert europepmc._read_zip(raw, "word/document.xml") == expected


SUPPLEMENT = _zip(
    {
        "S1.txt": b"2n = 4x\n=  40",
        "S2.html": b"<p>Chr <b>1A</b></p>",
        "S3.xlsx": _zip(
            {"xl/sharedStrings.xml": b'<sst><si><t>ploidy</t></si><si><t xml:space="p">4x</t></si></sst>'}
        ),
        "S4.docx": _zip({"word/document.xml": b"<w:p><w:t>tetraploid</w:t><w:t>line</w:t></w:p>"}),
        "S5.pdf": b"%PDF",
        "S6.csv": b"   ",
    }
)


@pytest.mark.parametrize(
    ("side_effect", "max_chars", "expected"),
    [
        pytest.param(
            [_response(content=SUPPLEMENT)],
            1000,
            "[Supplementary file: S1.txt]\n2n = 4x = 40\n\n"
            "[Supplementary file: S2.html]\nChr 1A\n\n"
            "[Supplementary file: S3.xlsx]\nploidy 4x\n\n"
            "[Supplementary file: S4.docx]\ntetraploid line",
            id="Readable files, pdf and blank skipped",
        ),
        pytest.param(
            [_response(content=SUPPLEMENT)], 5, "[Supplementary file: S1.txt]\n2n = ", id="Truncated"
        ),
        pytest.param([_response(content=_zip({"S1.pdf": b"%PDF"}))], 1000, None, id="Nothing readable"),
        pytest.param([_response(content=b"not a zip")], 1000, None, id="Bad zip"),
        pytest.param([_response(status=404, content=SUPPLEMENT)], 1000, None, id="Not found"),
        pytest.param(requests.ConnectionError("down"), 1000, None, id="Request error"),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.requests.get")
def test_get_supplementary_text_by_pmcid(
    mock_get: MagicMock, side_effect: Any, max_chars: int, expected: str | None
) -> None:
    """Test :func:`europepmc.get_supplementary_text_by_pmcid`.

    Args:
        mock_get: Patched `requests.get()`.
        side_effect: Response returned, or exception raised, by the request.
        max_chars: Maximum number of characters to keep.
        expected: Expected concatenated text.

    """
    mock_get.side_effect = side_effect
    assert europepmc.get_supplementary_text_by_pmcid("PMC1", max_chars=max_chars) == expected


@pytest.mark.parametrize(
    ("record", "expected"),
    [
        pytest.param({"abstractText": "Diploid."}, "Diploid.", id="Found"),
        pytest.param(None, None, id="Not found"),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.search_by_pmid")
def test_get_abstract_by_pmid(mock_pmid: MagicMock, record: dict | None, expected: str | None) -> None:
    """Test :func:`europepmc.get_abstract_by_pmid`.

    Args:
        mock_pmid: Patched `search_by_pmid()`.
        record: Record returned by the PMID lookup.
        expected: Expected abstract.

    """
    mock_pmid.return_value = record
    assert europepmc.get_abstract_by_pmid("1") == expected


@pytest.mark.parametrize(
    ("paper", "fulltext", "supplementary", "abstract", "expected"),
    [
        pytest.param(
            {"pmcid": "PMC1", "pmid": "1"},
            "<article/>",
            "[Supplementary file: S1.txt]\n2n",
            None,
            {
                "source": "fulltext",
                "pmcid": "PMC1",
                "text": "<article/>",
                "supplementary": "[Supplementary file: S1.txt]\n2n",
            },
            id="Full text with supplementary",
        ),
        pytest.param(
            {"pmcid": "PMC1"},
            "<article/>",
            None,
            None,
            {"source": "fulltext", "pmcid": "PMC1", "text": "<article/>"},
            id="Full text only",
        ),
        pytest.param(
            {"pmcid": "PMC1", "pmid": "1"},
            None,
            None,
            "Diploid.",
            {"source": "abstract", "pmid": "1", "text": "Diploid."},
            id="No full text, abstract by PMID",
        ),
        pytest.param(
            {"pmid": "1", "abstractText": "From search."},
            None,
            None,
            None,
            {"source": "search_result", "text": "From search."},
            id="Search-result abstract",
        ),
        pytest.param({}, None, None, None, {"source": "none", "text": None}, id="Nothing"),
    ],
)
@patch("ensembl.io.genomio.literature.europepmc.get_abstract_by_pmid")
@patch("ensembl.io.genomio.literature.europepmc.get_supplementary_text_by_pmcid")
@patch("ensembl.io.genomio.literature.europepmc.get_fulltext_by_pmcid")
def test_get_best_available_text(
    mock_fulltext: MagicMock,
    mock_supp: MagicMock,
    mock_abstract: MagicMock,
    paper: dict,
    fulltext: str | None,
    supplementary: str | None,
    abstract: str | None,
    expected: dict,
) -> None:
    """Test :func:`europepmc.get_best_available_text`.

    Args:
        mock_fulltext: Patched `get_fulltext_by_pmcid()`.
        mock_supp: Patched `get_supplementary_text_by_pmcid()`.
        mock_abstract: Patched `get_abstract_by_pmid()`.
        paper: Paper record.
        fulltext: Full text returned by PMCID.
        supplementary: Supplementary text returned by PMCID.
        abstract: Abstract returned by PMID.
        expected: Expected text record.

    """
    mock_fulltext.return_value = fulltext
    mock_supp.return_value = supplementary
    mock_abstract.return_value = abstract
    assert europepmc.get_best_available_text(paper) == expected
