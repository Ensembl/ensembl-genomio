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
"""Unit testing of `ensembl.io.genomio.literature.metadata` module.

Typical usage example::
    $ pytest test_metadata.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import requests

from ensembl.io.genomio.literature import metadata

REPORT = {
    "accession": "GCA_003086295.3",
    "assembly_info": {"assembly_name": "arahy.Tifrunner.gnm2", "biosample": {"linked_pmid": "31043759"}},
    "organism": {"tax_id": 3818, "organism_name": "Arachis hypogaea", "common_name": "peanut"},
    "assembly_stats": {"total_number_of_chromosomes": 20},
}


def _response(payload: Any = None, error: Exception | None = None) -> MagicMock:
    """Build a fake `requests.Response`.

    Args:
        payload: Value returned by ``json()``.
        error: Exception raised by ``json()``, if any.

    """
    response = MagicMock()
    response.json.side_effect = error
    response.json.return_value = payload
    return response


def _goat(attributes: dict) -> list[MagicMock]:
    """Wrap GoaT record attributes in a one-item list of fake responses.

    Args:
        attributes: The ``records[0].record.attributes`` mapping.

    """
    return [_response({"records": [{"record": {"attributes": attributes}}]})]


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param([_response({"reports": [REPORT]})], [REPORT], id="Reports"),
        pytest.param([_response({})], [], id="No reports key"),
        pytest.param(requests.ConnectionError("down"), [], id="Request error"),
        pytest.param([_response(error=ValueError("bad json"))], [], id="Invalid JSON"),
    ],
)
@patch("ensembl.io.genomio.literature.metadata.requests.get")
def test_fetch_assembly_report(mock_get: MagicMock, side_effect: Any, expected: list) -> None:
    """Test :func:`metadata._fetch_assembly_report`.

    Args:
        mock_get: Patched `requests.get()`.
        side_effect: What the patched call returns or raises.
        expected: Expected list of reports.

    """
    mock_get.side_effect = side_effect
    assert metadata._fetch_assembly_report("GCA_003086295.3") == expected


@pytest.mark.parametrize(
    ("accession", "reports", "expected"),
    [
        pytest.param(
            "GCA_003086295.3",
            [[REPORT]],
            {
                "assembly_accession": "GCA_003086295.3",
                "assembly_name": "arahy.Tifrunner.gnm2",
                "taxon_id": "3818",
                "scientific_name": "Arachis hypogaea",
                "common_name": "peanut",
                "linked_pmids": ["31043759"],
                "chromosome_number": 20,
                "chromosome_source": "ncbi",
            },
            id="Report found",
        ),
        pytest.param(
            "GCA_000001735.1",
            [
                [],
                [
                    {
                        "accession": "GCA_000001735.2",
                        "assembly_info": {"linked_pmid": ["1"], "chromosome_count": 5},
                    }
                ],
            ],
            {
                "assembly_accession": "GCA_000001735.2",
                "assembly_name": None,
                "taxon_id": "",
                "scientific_name": None,
                "common_name": None,
                "linked_pmids": ["1"],
                "chromosome_number": 5,
                "chromosome_source": "ncbi",
            },
            id="Stale version retried",
        ),
        pytest.param(
            "GCA_000001735",
            [[]],
            {"assembly_accession": "GCA_000001735", "chromosome_number": None, "chromosome_source": None},
            id="No report, no version to drop",
        ),
        pytest.param(
            "GCA_000001735.1",
            [[], []],
            {"assembly_accession": "GCA_000001735.1", "chromosome_number": None, "chromosome_source": None},
            id="No report after retry",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.metadata._fetch_assembly_report")
def test_fetch_assembly_metadata(
    mock_report: MagicMock, accession: str, reports: list, expected: dict
) -> None:
    """Test :func:`metadata.fetch_assembly_metadata`.

    Args:
        mock_report: Patched `_fetch_assembly_report()`.
        accession: Accession to look up.
        reports: Reports returned by each successive lookup.
        expected: Expected metadata dict.

    """
    mock_report.side_effect = reports
    assert metadata.fetch_assembly_metadata(accession) == expected


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        pytest.param(
            "<TaxaSet><Taxon><OtherNames><CommonName>groundnut</CommonName>"
            "<GenbankCommonName> peanut </GenbankCommonName></OtherNames></Taxon></TaxaSet>",
            "peanut",
            id="GenBank common name preferred",
        ),
        pytest.param(
            "<TaxaSet><Taxon><OtherNames><CommonName>groundnut</CommonName></OtherNames></Taxon></TaxaSet>",
            "groundnut",
            id="Common name",
        ),
        pytest.param("<TaxaSet><Taxon><OtherNames/></Taxon></TaxaSet>", None, id="No common name"),
        pytest.param("<TaxaSet/>", None, id="No taxon"),
        pytest.param("<TaxaSet>", None, id="Malformed XML"),
    ],
)
@patch("ensembl.io.genomio.literature.metadata.entrez_get")
def test_fetch_taxonomy_common_name(mock_entrez: MagicMock, xml: str, expected: str | None) -> None:
    """Test :func:`metadata.fetch_taxonomy_common_name`.

    Args:
        mock_entrez: Patched `entrez_get()`.
        xml: Taxonomy record returned by Entrez.
        expected: Expected common name.

    """
    mock_entrez.return_value.text = xml
    assert metadata.fetch_taxonomy_common_name("3818") == expected


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param(
            _goat(
                {
                    "ploidy": {"aggregation_source": "direct", "min": 4, "max": 4, "value": 4},
                    "chromosome_number": {"aggregation_source": "direct", "value": 40},
                }
            ),
            {"reference_ploidy": 4, "reference_chromosome": 40},
            id="Direct values",
        ),
        pytest.param(
            _goat(
                {
                    "ploidy": {"aggregation_source": "ancestor", "value": 4},
                    "ploidy_inferred": {"aggregation_source": "direct", "value": 2},
                    "chromosome_number": {"aggregation_source": "direct", "median": 10},
                }
            ),
            {"reference_ploidy": 2, "reference_chromosome": 10},
            id="Inferred value skipped, next key and value field used",
        ),
        pytest.param(
            _goat({"ploidy": {"aggregation_source": "direct", "min": 2, "max": 4}}),
            {"reference_ploidy": None, "reference_chromosome": None},
            id="Direct records disagree",
        ),
        pytest.param(
            [_response({"records": []})],
            {"reference_ploidy": None, "reference_chromosome": None},
            id="No records",
        ),
        pytest.param(
            requests.ConnectionError("down"),
            {"reference_ploidy": None, "reference_chromosome": None},
            id="Request error",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.metadata.requests.get")
def test_fetch_goat_reference(mock_get: MagicMock, side_effect: Any, expected: dict) -> None:
    """Test :func:`metadata.fetch_goat_reference`.

    Args:
        mock_get: Patched `requests.get()`.
        side_effect: Fake GoaT responses, or the exception the request raises.
        expected: Expected reference values.

    """
    mock_get.side_effect = side_effect
    assert metadata.fetch_goat_reference("3818") == expected


@pytest.mark.parametrize(
    ("assembly", "common_name", "expected", "looked_up"),
    [
        pytest.param({"taxon_id": ""}, None, {"taxon_id": ""}, False, id="No taxon ID"),
        pytest.param(
            {"taxon_id": "3818", "common_name": "peanut"},
            None,
            {"taxon_id": "3818", "common_name": "peanut", "reference_ploidy": 4, "reference_chromosome": 20},
            False,
            id="Common name already set",
        ),
        pytest.param(
            {"taxon_id": "3818"},
            "peanut",
            {
                "taxon_id": "3818",
                "common_name": "peanut",
                "common_name_source": "ncbi_taxonomy",
                "reference_ploidy": 4,
                "reference_chromosome": 20,
            },
            True,
            id="Common name filled from Taxonomy",
        ),
        pytest.param(
            {"taxon_id": "3818"},
            None,
            {"taxon_id": "3818", "reference_ploidy": 4, "reference_chromosome": 20},
            True,
            id="Common name not found",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.metadata.fetch_goat_reference")
@patch("ensembl.io.genomio.literature.metadata.fetch_taxonomy_common_name")
def test_enrich_assembly_metadata(
    mock_common: MagicMock,
    mock_goat: MagicMock,
    assembly: dict,
    common_name: str | None,
    expected: dict,
    looked_up: bool,
) -> None:
    """Test :func:`metadata.enrich_assembly_metadata`.

    Args:
        mock_common: Patched `fetch_taxonomy_common_name()`.
        mock_goat: Patched `fetch_goat_reference()`.
        assembly: Assembly metadata to enrich.
        common_name: Common name returned by NCBI Taxonomy.
        expected: Expected enriched metadata.
        looked_up: Whether NCBI Taxonomy is expected to be queried.

    """
    mock_common.return_value = common_name
    mock_goat.return_value = {"reference_ploidy": 4, "reference_chromosome": 20}
    assert metadata.enrich_assembly_metadata(assembly) == expected
    assert mock_common.called == looked_up
