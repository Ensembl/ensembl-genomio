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
"""Unit testing of `ensembl.io.genomio.literature.parse` module.

Typical usage example::
    $ pytest test_parse.py

"""

from xml.etree import ElementTree as ET

import pytest

from ensembl.io.genomio.literature import parse

ARTICLE_XML = """<article>
<front><abstract><p>Wheat is <italic>hexaploid</italic>.</p></abstract></front>
<body>
<sec><title>Methods</title><p>A longer methods paragraph.</p></sec>
<sec><title>Methods</title><p>Short.</p></sec>
<sec><title>Results</title><p>2n = 6x = 42</p></sec>
<sec><p>No title, skipped.</p></sec>
</body>
<table-wrap><label>Table 1</label><caption><p>Karyotype</p></caption>
<table><thead><tr><th>Chr</th><th>Length</th></tr></thead>
<tbody><tr><td>1A</td><td>594</td></tr><tr><td></td><td></td></tr></tbody></table></table-wrap>
</article>"""


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param("", "", id="Empty"),
        pytest.param("2n = 4x  =\n30", "2n = 4x = 30", id="Unicode whitespace"),
        pytest.param(r"ploidy $4x$ \cite{ref} here", "ploidy 4x here", id="TeX commands and math"),
        pytest.param(r"a \documentclass{x} body \end{document} b", "a b", id="TeX document block"),
    ],
)
def test_clean_text(text: str, expected: str) -> None:
    """Test :func:`parse.clean_text`.

    Args:
        text: Raw input text.
        expected: Expected cleaned text.

    """
    assert parse.clean_text(text) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param("  <article></article>", "xml", id="XML"),
        pytest.param("Plain abstract.", "abstract", id="Abstract"),
    ],
)
def test_detect_text_type(text: str, expected: str) -> None:
    """Test :func:`parse.detect_text_type`.

    Args:
        text: Raw input text.
        expected: Expected text type.

    """
    assert parse.detect_text_type(text) == expected


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        pytest.param(
            "<a><table-wrap><label>Table 1</label><caption><p>Karyotype</p></caption><table>"
            "<thead><tr><th>Chr</th><th>Length</th></tr></thead>"
            "<tbody><tr><td>1A</td><td>594</td></tr><tr><td></td><td></td></tr></tbody>"
            "</table></table-wrap></a>",
            {"Table 1": "Karyotype Chr: 1A | Length: 594"},
            id="Header and rows",
        ),
        pytest.param(
            "<a><table-wrap><table><tr><td>1</td><td>2</td><td>3</td></tr>"
            "<tr><td>x</td><td>y</td></tr></table></table-wrap></a>",
            {"Table 1": "1: 1 | 2: 2 | 3: 3 ; x | y"},
            id="No thead, row length differs from header",
        ),
        pytest.param(
            "<a><table-wrap><label>T</label><caption><p>Long caption</p></caption></table-wrap>"
            "<table-wrap><label>T</label><caption><p>Short</p></caption></table-wrap></a>",
            {"T": "Long caption"},
            id="Repeated title keeps longest body",
        ),
        pytest.param(
            "<a><table-wrap><table>"
            + "".join(f"<tr><td>{i}</td></tr>" for i in range(45))
            + "</table></table-wrap></a>",
            {"Table 1": " ; ".join(f"0: {i}" for i in range(40))},
            id="Only the first 40 rows",
        ),
        pytest.param("<a><table-wrap><label>Empty</label></table-wrap></a>", {}, id="Empty table skipped"),
    ],
)
def test_parse_xml_tables(xml: str, expected: dict[str, str]) -> None:
    """Test :func:`parse.parse_xml_tables`.

    Args:
        xml: XML document holding one or more ``table-wrap`` elements.
        expected: Expected ``{title: body}`` sections.

    """
    assert parse.parse_xml_tables(ET.fromstring(xml)) == expected


@pytest.mark.parametrize(
    ("xml_text", "expected"),
    [
        pytest.param(
            ARTICLE_XML,
            {
                "Abstract": "Wheat is hexaploid.",
                "Methods": "A longer methods paragraph.",
                "Results": "2n = 6x = 42",
                "Table 1": "Karyotype Chr: 1A | Length: 594",
            },
            id="Full article",
        ),
        pytest.param("<article><unclosed></article>", {}, id="Malformed XML"),
    ],
)
def test_parse_xml_sections(xml_text: str, expected: dict[str, str]) -> None:
    """Test :func:`parse.parse_xml_sections`.

    Args:
        xml_text: PMC full-text XML.
        expected: Expected ``{section title: text}`` dict.

    """
    assert parse.parse_xml_sections(xml_text) == expected


def test_parse_abstract() -> None:
    """Test :func:`parse.parse_abstract`."""
    assert parse.parse_abstract("<title>Abstract</title> <p>Diploid  2n = 2x</p>") == {
        "Abstract": "Abstract Diploid 2n = 2x"
    }
    assert parse.parse_abstract(None) == {"Abstract": ""}  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("supp_text", "expected"),
    [
        pytest.param("", {}, id="Empty"),
        pytest.param(
            "[Supplementary file: s1.txt] Table S1 [Supplementary file: s2.txt]  "
            "[Supplementary file: s3.txt] x",
            {"Supplementary: s1.txt": "Table S1", "Supplementary: s3.txt": "x"},
            id="Markers, empty block skipped",
        ),
        pytest.param("No markers here", {"Supplementary": "No markers here"}, id="No markers"),
    ],
)
def test_parse_supplementary(supp_text: str, expected: dict[str, str]) -> None:
    """Test :func:`parse.parse_supplementary`.

    Args:
        supp_text: Concatenated supplementary text.
        expected: Expected sections.

    """
    assert parse.parse_supplementary(supp_text) == expected


def test_split_into_chunks() -> None:
    """Test :func:`parse.split_into_chunks`."""
    sections = {"A": "One two. Three four. Five.", "B": "Short."}
    assert parse.split_into_chunks(sections, chunk_size=12) == [
        {"section": "A", "text": "One two.", "chunk_index": 0},
        {"section": "A", "text": "Three four.", "chunk_index": 1},
        {"section": "A", "text": "Five.", "chunk_index": 2},
        {"section": "B", "text": "Short.", "chunk_index": 0},
    ]


@pytest.mark.parametrize(
    ("text_data", "expected_sections"),
    [
        pytest.param({}, {}, id="No text"),
        pytest.param(
            {"source": "pmc", "text": "<article><sec><title>M</title><p>Body.</p></sec></article>"},
            {"M": "Body."},
            id="Full-text XML",
        ),
        pytest.param(
            {"source": "epmc", "text": "<title>Abstract</title> <p>Tetraploid.</p>"},
            {"Abstract": "Abstract Tetraploid."},
            id="Tagged abstract falls back",
        ),
        pytest.param(
            {"source": "epmc", "text": "Diploid.", "supplementary": "[Supplementary file: s1] Chr 1A"},
            {"Abstract": "Diploid.", "Supplementary: s1": "Chr 1A"},
            id="Abstract with supplementary",
        ),
    ],
)
def test_parse_paper(text_data: dict, expected_sections: dict[str, str]) -> None:
    """Test :func:`parse.parse_paper`.

    Args:
        text_data: Paper text record as produced by the fetch step.
        expected_sections: Expected parsed sections.

    """
    result = parse.parse_paper(text_data)
    assert result["source"] == text_data.get("source", "none")
    assert result["sections"] == expected_sections
    assert result["chunks"] == parse.split_into_chunks(expected_sections)
    assert result["n_sections"] == len(expected_sections)
    assert result["n_chunks"] == len(result["chunks"])
