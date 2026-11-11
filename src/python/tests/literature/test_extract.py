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
"""Unit testing of `ensembl.io.genomio.literature.extract` module.

Typical usage example::
    $ pytest test_extract.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

from unittest.mock import MagicMock, patch

import pytest

from ensembl.io.genomio.literature import extract

PEANUT_NAMES = ["arachis", "hypogaea", "peanut"]


@pytest.mark.parametrize(
    ("section", "expected"),
    [
        pytest.param("Title", 2.5, id="Title"),
        pytest.param("Materials and Methods", 0.8, id="Methods"),
        pytest.param("Data Availability Statement", 0.0, id="Noise section"),
        pytest.param("Karyotype", 1.0, id="Unmatched"),
    ],
)
def test_get_section_weight(section: str, expected: float) -> None:
    """Test :func:`extract.get_section_weight`.

    Args:
        section: Section title.
        expected: Expected weight.

    """
    assert extract.get_section_weight(section) == expected


@pytest.mark.parametrize(
    ("scientific_name", "common_name", "expected"),
    [
        pytest.param(
            "Triticum aestivum", "bread wheat", ["triticum", "aestivum", "wheat"], id="Head noun only"
        ),
        pytest.param("Oryza sativa var. japonica", None, ["oryza", "sativa"], id="Binomial of a variety"),
        pytest.param("Hippophae rhamnoides", "sea-buckthorn", ["hippophae", "rhamnoides", "buckthorn"]),
        pytest.param("Bos taurus", "ox", ["bos", "taurus"], id="Short common name"),
        pytest.param(None, None, [], id="Nothing"),
    ],
)
def test_target_names(scientific_name: str | None, common_name: str | None, expected: list[str]) -> None:
    """Test :func:`extract.target_names`.

    Args:
        scientific_name: Scientific name.
        common_name: Common name.
        expected: Expected lower-cased names.

    """
    assert extract.target_names(scientific_name, common_name) == expected


@pytest.mark.parametrize(
    ("context", "names", "expected"),
    [
        pytest.param("Dioscorea nipponica is tetraploid", [], False, id="No names to check"),
        pytest.param("Unlike D. nipponica, peanut is tetraploid", PEANUT_NAMES, False, id="Target named"),
        pytest.param(
            "as in D. nipponica, which is tetraploid", PEANUT_NAMES, True, id="Abbreviated other species"
        ),
        pytest.param("A. duranensis is diploid", PEANUT_NAMES, True, id="Same initial, other epithet"),
        pytest.param("Dioscorea nipponica is tetraploid", PEANUT_NAMES, True, id="Full binomial"),
        pytest.param("Our data show a tetraploid", PEANUT_NAMES, False, id="Prose, not a species"),
    ],
)
def test_names_other_species(context: str, names: list[str], expected: bool) -> None:
    """Test :func:`extract._names_other_species`.

    Args:
        context: Text around a ploidy mention.
        names: Target organism names.
        expected: Whether the context names only another species.

    """
    assert extract._names_other_species(context, names) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param("Yeast has haploid and diploid phases; haploids mate.", 2, id="Organism-level"),
        pytest.param("A doubled haploid line and a pseudo-haploid.", 0, id="Doubled and pseudo haploid"),
        pytest.param("The haploid genome size and haploid chromosome number.", 0, id="Technical"),
        pytest.param("A dihaploid potato.", 0, id="Dihaploid"),
    ],
)
def test_count_haploid_mentions(text: str, expected: int) -> None:
    """Test :func:`extract.count_haploid_mentions`.

    Args:
        text: Text to scan.
        expected: Expected number of organism-level mentions.

    """
    assert extract.count_haploid_mentions(text) == expected


@pytest.mark.parametrize(
    ("text", "names", "expected"),
    [
        pytest.param("x is tetraploid", [], True, id="No names to check"),
        pytest.param("Here, our tetraploid genome is presented.", PEANUT_NAMES, True, id="Self-reference"),
        pytest.param(
            "Peanut " + "x" * 100 + " is tetraploid", PEANUT_NAMES, True, id="Name within the window"
        ),
        pytest.param("Peanut " + "x" * 200 + " is tetraploid", PEANUT_NAMES, False, id="Name too far"),
    ],
)
def test_is_anchored(text: str, names: list[str], expected: bool) -> None:
    """Test :func:`extract._is_anchored` on the mention of "tetraploid".

    Args:
        text: Text holding the mention.
        names: Target organism names.
        expected: Whether the mention is anchored to the target.

    """
    start = text.index("tetraploid")
    assert extract._is_anchored(text, start, start + len("tetraploid"), names) == expected


@pytest.mark.parametrize(
    ("level", "two_n", "count", "expected"),
    [
        pytest.param(4, 40, None, True, id="No count to check"),
        pytest.param(4, 40, 40, True, id="Full set"),
        pytest.param(4, 40, 20, True, id="Haploid set"),
        pytest.param(4, 40, 22, True, id="Haploid set plus organelles"),
        pytest.param(4, 38, 19, True, id="Haploid set, aneuploid formula"),
        pytest.param(6, 42, 14, True, id="Multiple of the base number"),
        pytest.param(6, 42, 49, False, id="Multiple above 2n"),
        pytest.param(4, 40, 5, False, id="Another organism"),
        pytest.param(3, 40, 26, False, id="Non-integer base number"),
    ],
)
def test_formula_matches_chromosomes(level: int, two_n: int, count: int | None, expected: bool) -> None:
    """Test :func:`extract.formula_matches_chromosomes`.

    Args:
        level: Ploidy level from the formula.
        two_n: Chromosome count from the formula.
        count: Assembled chromosome count.
        expected: Whether the formula fits the assembly.

    """
    assert extract.formula_matches_chromosomes(level, two_n, count) == expected


@pytest.mark.parametrize(
    ("text", "count", "names", "expected"),
    [
        pytest.param("Peanut has 2n = 4x = 40.", 20, [], [(4, True, False, False, True, True)], id="Formula"),
        pytest.param("Yam has 2n = 4x = 40.", 5, [], [], id="Formula inconsistent with the assembly"),
        pytest.param("2n = 30x = 60", None, [], [], id="Formula level out of range"),
        pytest.param(
            "An allotetraploid derived from two diploid species.",
            None,
            [],
            [(4, True, False, False, True, False), (2, False, True, False, True, False)],
            id="Compound term, ancestor context",
        ),
        pytest.param(
            "A diploid-aware assembler.", None, [], [(2, False, False, True, True, False)], id="Technical"
        ),
        pytest.param(
            "Dioscorea nipponica is tetraploid.",
            None,
            PEANUT_NAMES,
            [(4, False, False, True, False, False)],
            id="Other species",
        ),
        pytest.param("Peanut is tetraploid.", None, PEANUT_NAMES, [(4, False, False, False, True, False)]),
        pytest.param("An amphidiploid; diploidization.", None, [], [], id="Not level words"),
    ],
)
def test_find_level_mentions(text: str, count: int | None, names: list[str], expected: list[tuple]) -> None:
    """Test :func:`extract.find_level_mentions` classification.

    Args:
        text: Text to scan.
        count: Assembled chromosome count.
        names: Target organism names.
        expected: Expected ``(level, is_compound, is_ancestor, is_off_target, is_anchored, is_formula)``
            for each mention.

    """
    mentions = extract.find_level_mentions(text, chromosome_count=count, names=names)
    fields = ("level", "is_compound", "is_ancestor", "is_off_target", "is_anchored", "is_formula")
    assert [tuple(mention[field] for field in fields) for mention in mentions] == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param("An allotetraploid, or amphidiploid, crop.", "amphidiploid", id="Most specific wins"),
        pytest.param("A segmental allopolyploid.", "segmental_allopolyploid", id="Segmental"),
        pytest.param("An auto-tetraploid potato.", "autopolyploid", id="Hyphenated"),
        pytest.param("Its allopolyploid ancestor.", None, id="Qualifies an ancestor"),
        pytest.param("A diploid.", None, id="None"),
        pytest.param(
            "Peanut arose from an allopolyploid ancestor. Today, peanut is an allotetraploid crop.",
            "allopolyploid",
            id="Later organism-level mention",
        ),
    ],
)
def test_detect_mechanism(text: str, expected: str | None) -> None:
    """Test :func:`extract.detect_mechanism`.

    Args:
        text: Text to scan.
        expected: Expected mechanism.

    """
    assert extract.detect_mechanism(text) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param("An aneuploid line.", "aneuploid", id="Aneuploid"),
        pytest.param("Endoreduplication in leaves.", "endopolyploid", id="Endopolyploid"),
        pytest.param("A diploid.", False, id="None"),
    ],
)
def test_detect_irregular(text: str, expected: str | bool) -> None:
    """Test :func:`extract.detect_irregular`.

    Args:
        text: Text to scan.
        expected: Expected irregularity term.

    """
    assert extract.detect_irregular(text) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param("A doubled-haploid line.", "doubled_haploid", id="Doubled haploid"),
        pytest.param("A double monoploid line.", "doubled_monoploid", id="Doubled monoploid"),
        pytest.param("A colchicine-induced line.", "induced_polyploid", id="Induced"),
        pytest.param("A diploid.", None, id="None"),
    ],
)
def test_detect_derivation(text: str, expected: str | None) -> None:
    """Test :func:`extract.detect_derivation`.

    Args:
        text: Text to scan.
        expected: Expected derivation.

    """
    assert extract.detect_derivation(text) == expected


@pytest.mark.parametrize(
    ("segments", "expected"),
    [
        pytest.param([], (None, None, 0.0, "no_text_evidence", "not_stated_in_paper", 0), id="No evidence"),
        pytest.param(
            [("Peanut is tetraploid.", 1.0, False)],
            (4, None, 0.5, "text_weighted_vote", "found_in_paper", 1),
            id="One mention, capped",
        ),
        pytest.param(
            [("Peanut is allotetraploid.", 1.0, False)],
            (4, "allopolyploid", 0.7, "text_weighted_vote", "found_in_paper", 1),
            id="Compound term counts twice",
        ),
        pytest.param(
            [("Peanut (2n = 4x = 40) is often called diploid.", 1.0, True)],
            (4, None, 0.95, "chromosome_formula", "found_in_paper", 2),
            id="Formula outweighs words, uncapped",
        ),
        pytest.param(
            [
                (
                    "Peanut is tetraploid. It is grown in warm regions across the world. "
                    "Its diploid progenitor is wild.",
                    1.0,
                    False,
                )
            ],
            (4, None, 0.5, "text_weighted_vote", "found_in_paper", 1),
            id="Ancestor mention not counted",
        ),
        pytest.param(
            [("An allotetraploid crop. Diploid. Diploid.", 1.0, False)],
            (4, "allopolyploid", 0.7, "text_weighted_vote", "found_in_paper", 3),
            id="Compound term outweighs plain words",
        ),
        pytest.param(
            [("A tetraploid.", 2.0, True), ("Diploid. Diploid. Diploid.", 0.8, False)],
            (2, None, 0.55, "text_weighted_vote", "found_in_paper", 4),
            id="Weighted vote",
        ),
        pytest.param(
            [("A tetraploid.", 0.0, True), ("", 1.0, False), ("A diploid.", 1.0, False)],
            (2, None, 0.5, "text_weighted_vote", "found_in_paper", 1),
            id="Noise and empty segments skipped",
        ),
        pytest.param(
            [("An allopolyploid species.", 1.0, False)],
            (None, "allopolyploid", 0.5, "no_text_evidence", "found_in_paper", 0),
            id="Mechanism without a level",
        ),
        pytest.param(
            [("This diploid line is not autopolyploid.", 1.0, False)],
            (2, None, 0.5, "text_weighted_vote", "found_in_paper", 1),
            id="Diploid has no mechanism",
        ),
    ],
)
def test_resolve_ploidy_fields(segments: list[tuple], expected: tuple) -> None:
    """Test :func:`extract.resolve_ploidy_fields` level vote, confidence and schema normalisation.

    Args:
        segments: Weighted ``(text, weight, is_priority)`` segments.
        expected: Expected ``(level, mechanism, confidence, method, status, evidence_count)``.

    """
    result = extract.resolve_ploidy_fields(segments)
    fields = ("level", "mechanism", "confidence", "method", "status", "evidence_count")
    assert tuple(result[field] for field in fields) == expected


@pytest.mark.parametrize(
    ("segments", "expected"),
    [
        pytest.param(
            [
                ("A doubled haploid line of a diploid, aneuploid crop.", 2.0, True),
                ("Haploid spores germinate. The line is diploid, not mixoploid.", 1.0, False),
                ("Zea mays is tetraploid in other work.", 1.0, False),
            ],
            ("aneuploid", "doubled_haploid", "doubled_haploid", 1, False, {2: 3.0, 4: 1.0}),
            id="Derivation in a priority section, winner not anchored",
        ),
        pytest.param(
            [("A diploid crop.", 2.0, True), ("Zea mays, a doubled haploid line, is diploid.", 1.0, False)],
            (False, "doubled_haploid", None, 0, True, {2: 3.0}),
            id="Derivation outside priority sections, winner anchored",
        ),
    ],
)
def test_resolve_ploidy_fields_guard_signals(segments: list[tuple], expected: tuple) -> None:
    """Test :func:`extract.resolve_ploidy_fields` reports the signals for the pipeline's precision guard.

    Args:
        segments: Weighted ``(text, weight, is_priority)`` segments about a diploid maize line.
        expected: Expected ``(irregular, derivation, priority_derivation, haploid_mentions, anchored,
            weighted_counts)``.

    """
    result = extract.resolve_ploidy_fields(segments, chromosome_count=20, names=["zea", "mays"])
    assert result["level"] == 2
    fields = (
        "irregular",
        "derivation",
        "priority_derivation",
        "haploid_mentions",
        "anchored",
        "weighted_counts",
    )
    assert tuple(result[field] for field in fields) == expected


@patch("ensembl.io.genomio.literature.extract.resolve_ploidy_fields")
def test_extract_ploidy(mock_resolve: MagicMock) -> None:
    """Test :func:`extract.extract_ploidy` marks priority sections and passes the checks through.

    Args:
        mock_resolve: Patched `resolve_ploidy_fields()`.

    """
    sections = {"Abstract": ("a", 2.0), "Methods": ("m", 0.8), "Background": ("b", 1.5)}
    assert extract.extract_ploidy(sections, chromosome_count=20, names=["zea"]) is mock_resolve.return_value
    mock_resolve.assert_called_once_with(
        [("a", 2.0, True), ("m", 0.8, False), ("b", 1.5, True)], chromosome_count=20, names=["zea"]
    )


@pytest.mark.parametrize(
    ("chromosome_number", "expected"),
    [
        pytest.param({"value": "20", "source": "ncbi"}, 20, id="From NCBI"),
        pytest.param({"value": "20", "source": "text_extraction"}, None, id="From text"),
        pytest.param({"value": "twenty", "source": "ncbi"}, None, id="Not a number"),
        pytest.param(None, None, id="Missing"),
    ],
)
def test_ncbi_chromosome_count(chromosome_number: dict | None, expected: int | None) -> None:
    """Test :func:`extract.ncbi_chromosome_count`.

    Args:
        chromosome_number: Chromosome-number record.
        expected: Expected count.

    """
    assert extract.ncbi_chromosome_count(chromosome_number) == expected


@pytest.mark.parametrize(
    ("sections", "expected"),
    [
        pytest.param(
            {"Methods": ("It has 2n = 4x = 52 chromosomes.", 0.8), "Abstract": ("We found 2n = 48.", 2.0)},
            ("48", "Abstract"),
            id="Highest-weight section first",
        ),
        pytest.param({"Results": ("2n = 4x = 52", 1.2)}, ("52", "Results"), id="Count, not the ploidy"),
        pytest.param({"Results": ("It has 24 chromosomes.", 1.2)}, ("24", "Results"), id="N chromosomes"),
        pytest.param(
            {"Results": ("a chromosome number of 7", 1.2)}, ("7", "Results"), id="Chromosome number of"
        ),
        pytest.param({"Data availability": ("2n = 48", 0.0)}, (None, None), id="Noise section skipped"),
        pytest.param(
            {"Methods": ("The k-mer size used for estimation = 21 bp.", 0.8)},
            (None, None),
            id="Equation in prose",
        ),
    ],
)
def test_extract_chromosome_number(sections: dict, expected: tuple) -> None:
    """Test :func:`extract.extract_chromosome_number`.

    Args:
        sections: Weighted sections.
        expected: Expected ``(value, section)``.

    """
    result = extract.extract_chromosome_number(sections)
    assert (result["value"], result.get("section")) == expected
    assert result["source"] == "text_extraction"


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param(
            "We used cultivar 'Tifrunner' and cv. 'Chinese Spring'.",
            ["Chinese Spring", "Tifrunner"],
            id="Quoted",
        ),
        pytest.param("Cultivar: Col-0 and strain TM-1.", ["Col-0", "TM-1"], id="Labelled"),
        pytest.param("Shown in cultivar Table and strain Fasta.", [], id="Stopwords"),
        pytest.param("Seeds of cultivar 'the reference line' were used.", [], id="Not a name"),
        pytest.param(
            "The cv. Nipponbare was sequenced.",
            ["Nipponbare"],
            id="Abbreviated cultivar mid-sentence",
        ),
    ],
)
def test_extract_cultivars(text: str, expected: list[str]) -> None:
    """Test :func:`extract.extract_cultivars`.

    Args:
        text: Paper text.
        expected: Expected sorted cultivar names.

    """
    assert extract.extract_cultivars(text) == expected


@pytest.mark.parametrize(
    ("sections", "expected"),
    [
        pytest.param(
            {"Methods": ("A female plant was used.", 0.8)}, ("female", "Methods"), id="Sequenced individual"
        ),
        pytest.param({"Results": ("Sex: male", 1.2)}, ("male", "Results"), id="Labelled"),
        pytest.param(
            {"Data availability": ("female individual", 0.0)}, ("unknown", None), id="Noise section"
        ),
        pytest.param(
            {"Methods": ("The maternal haplotype was assembled.", 0.8)},
            ("unknown", None),
            id="Haplotype parent is not a sex",
        ),
    ],
)
def test_extract_sex(sections: dict, expected: tuple) -> None:
    """Test :func:`extract.extract_sex`.

    Args:
        sections: Weighted sections.
        expected: Expected ``(value, section)``.

    """
    result = extract.extract_sex(sections)
    assert (result["value"], result.get("section")) == expected


@pytest.mark.parametrize(
    ("abstract", "sections", "expected"),
    [
        pytest.param(
            "Apple orchards of Malus domestica Borkh. were sampled.",
            {},
            "Malus domestica",
            id="Abstract, name with authority preferred",
        ),
        pytest.param(
            "",
            {"Methods": ("Prunus persica", 0.8), "Introduction": ("Pyrus communis grows", 1.5)},
            "Pyrus communis",
            id="High-weight section first",
        ),
        pytest.param("", {"Methods": ("Prunus persica", 0.8)}, "Prunus persica", id="Any non-noise section"),
        pytest.param("", {"Methods": ("2n = 40", 0.8)}, "unknown", id="Not found"),
        pytest.param(
            "The genome of Malus domestica is presented.",
            {},
            "Malus domestica",
            id="Binomial after prose",
        ),
    ],
)
def test_extract_species(abstract: str, sections: dict, expected: str) -> None:
    """Test :func:`extract.extract_species`.

    Args:
        abstract: Abstract text.
        sections: Weighted sections.
        expected: Expected species name.

    """
    assert extract.extract_species({"Abstract": abstract}, sections) == expected


@pytest.mark.parametrize(
    ("assembly", "expected_species", "expected_chromosomes", "expected_count"),
    [
        pytest.param(
            {
                "assembly_accession": "GCA_003086295.3",
                "assembly_name": "arahy",
                "taxon_id": "3818",
                "scientific_name": "Arachis hypogaea",
                "common_name": "peanut",
                "chromosome_number": 20,
                "chromosome_source": "ncbi",
            },
            {"value": "Arachis hypogaea", "source": "ncbi"},
            {"value": "20", "source": "ncbi"},
            20,
            id="NCBI fields used",
        ),
        pytest.param(
            {"assembly_accession": "GCA_1.1"},
            {"value": "Pyrus communis", "source": "text_extraction"},
            {"value": "38", "section": "Abstract", "source": "text_extraction"},
            None,
            id="Text fills the gaps",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.extract.extract_ploidy")
def test_extract_metadata(
    mock_ploidy: MagicMock,
    assembly: dict,
    expected_species: dict,
    expected_chromosomes: dict,
    expected_count: int | None,
) -> None:
    """Test :func:`extract.extract_metadata` prefers NCBI fields, adds the title and skips noise sections.

    Args:
        mock_ploidy: Patched `extract_ploidy()`, returning a paper-derived level.
        assembly: Assembly metadata.
        expected_species: Expected species record.
        expected_chromosomes: Expected chromosome-number record (subset).
        expected_count: Chromosome count expected to be passed to the ploidy check.

    """
    mock_ploidy.return_value = {"level": 4, "tuple": [4, None, False, None], "confidence": 0.9, "method": "m"}
    parsed = {
        "sections": {
            "Abstract": "Pyrus communis with 2n = 38 chromosomes and a female individual.",
            "Methods": "We used cultivar 'Bartlett'.",
            "Data availability": "cultivar 'Hidden'",
        }
    }
    result = extract.extract_metadata(parsed, assembly, paper_title="A tetraploid genome")
    (weighted,) = mock_ploidy.call_args.args
    assert weighted["Title"] == ("A tetraploid genome", 2.5)
    assert mock_ploidy.call_args.kwargs["chromosome_count"] == expected_count
    assert "Title" not in parsed["sections"]
    assert result["species"] == expected_species
    assert expected_chromosomes.items() <= result["chromosome_number"].items()
    assert (result["cultivars"], result["sex"]["value"]) == (["Bartlett"], "female")


@pytest.mark.parametrize(
    ("assembly", "expected"),
    [
        pytest.param(
            {"taxon_id": "3818", "reference_ploidy": 6}, (6, "GoaT", "from_reference"), id="GoaT preferred"
        ),
        pytest.param({"taxon_id": "3818"}, (4, "reference_db", "from_reference"), id="Built-in table"),
        pytest.param({"taxon_id": "9606"}, (None, None, None), id="No reference"),
    ],
)
@patch("ensembl.io.genomio.literature.extract.extract_ploidy")
def test_extract_metadata_reference_fallback(mock_ploidy: MagicMock, assembly: dict, expected: tuple) -> None:
    """Test :func:`extract.extract_metadata` uses a sourced reference ploidy when the paper is silent.

    Args:
        mock_ploidy: Patched `extract_ploidy()`, returning no level.
        assembly: Assembly metadata.
        expected: Expected ``(level, source, status)``.

    """
    mock_ploidy.return_value = {
        "level": None,
        "tuple": [None, None, False, None],
        "confidence": 0.0,
        "method": "m",
    }
    ploidy = extract.extract_metadata({"sections": {"Abstract": "No ploidy here."}}, assembly)["ploidy"]
    assert (ploidy["level"], ploidy.get("source"), ploidy.get("status")) == expected
