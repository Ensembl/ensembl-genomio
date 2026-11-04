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
"""Unit testing of `ensembl.io.genomio.literature.pipeline` module.

Typical usage example::
    $ pytest test_pipeline.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

from contextlib import nullcontext as does_not_raise
import json
import os
from pathlib import Path
from typing import Any, ContextManager
from unittest.mock import MagicMock, call, patch

import pytest

from ensembl.io.genomio import literature
from ensembl.io.genomio.literature import pipeline

ACCESSION = "GCA_003086295.3"
ASSEMBLY = {
    "assembly_accession": ACCESSION,
    "assembly_name": "arahy",
    "taxon_id": "3818",
    "scientific_name": "Arachis hypogaea",
    "chromosome_number": 20,
    "chromosome_source": "ncbi",
    "reference_ploidy": 4,
}


def _paper(retrieval_source: str = "elink", text: str = "Peanut is tetraploid.") -> dict:
    """Build a prepared paper record as returned by `fetch_papers_for_assembly()`.

    Args:
        retrieval_source: How the paper was found.
        text: Paper text.

    """
    return {
        "title": "Peanut genome",
        "pmcid": "PMC1",
        "pmid": "1",
        "doi": "10.1/x",
        "retrieval_source": retrieval_source,
        "relevance_score": 9.0,
        "validation": {"pmcid_valid": True, "pmid_valid": True},
        "text_data": {"source": "fulltext", "text": text},
    }


@pytest.mark.parametrize(
    ("paper", "expected"),
    [
        pytest.param(
            {"title": "Genome of peanut", "text_data": {"text": "Nuclear genome."}}, True, id="Nuclear"
        ),
        pytest.param({"title": "The Plastome of peanut"}, False, id="Organelle in title"),
        pytest.param(
            {"title": "Peanut", "text_data": {"text": "Its mtDNA was assembled."}}, False, id="In abstract"
        ),
        pytest.param(
            {"title": "Peanut", "text_data": {"text": "x" * 500 + " chloroplast"}}, True, id="Beyond 500"
        ),
    ],
)
def test_is_nuclear_genome_paper(paper: dict, expected: bool) -> None:
    """Test :func:`pipeline.is_nuclear_genome_paper`.

    Args:
        paper: Paper record.
        expected: Whether the paper is about the nuclear genome.

    """
    assert pipeline.is_nuclear_genome_paper(paper) == expected


@pytest.mark.parametrize(
    ("papers", "expected_id"),
    [
        pytest.param([], None, id="No papers"),
        pytest.param([{"id": "a", "text_data": {"source": "none"}}, {"id": "b"}], "a", id="No usable text"),
        pytest.param(
            [
                {"id": "a", "relevance_score": 3.0, "text_data": {"source": "fulltext"}},
                {"id": "b", "relevance_score": 5.0, "text_data": {"source": "search_result"}},
                {"id": "c", "relevance_score": 9.0, "text_data": {"source": "none"}},
            ],
            "b",
            id="Highest score with usable text",
        ),
        pytest.param(
            [
                {"id": "a", "relevance_score": 5.0, "text_data": {"source": "abstract"}},
                {"id": "b", "relevance_score": 5.0, "text_data": {"source": "fulltext"}},
            ],
            "b",
            id="Tie goes to full text",
        ),
    ],
)
def test_select_best_paper(papers: list[dict], expected_id: str | None) -> None:
    """Test :func:`pipeline.select_best_paper`.

    Args:
        papers: Candidate papers.
        expected_id: ID of the expected paper.

    """
    best = pipeline.select_best_paper(papers)
    assert (best or {}).get("id") == expected_id


@pytest.mark.parametrize(
    ("paper_level", "goat", "expected"),
    [
        pytest.param(4, None, "no_reference", id="No reference"),
        pytest.param(None, 4, "no_paper_ploidy", id="No paper ploidy"),
        pytest.param(4, 4, "agree", id="Agree"),
        pytest.param(2, 4, "disagree", id="Disagree"),
    ],
)
def test_goat_agreement(paper_level: int | None, goat: int | None, expected: str) -> None:
    """Test :func:`pipeline._goat_agreement`.

    Args:
        paper_level: Paper-derived ploidy level.
        goat: GoaT reference ploidy.
        expected: Expected agreement flag.

    """
    assert pipeline._goat_agreement(paper_level, goat) == {"goat_ploidy": goat, "agreement": expected}


SOUND = {
    "level": 4,
    "status": "found_in_paper",
    "method": "consensus",
    "confidence": 0.8,
    "evidence": "Peanut is tetraploid.",
    "anchored": True,
    "haploid_mentions": 0,
    "reference_check": {"goat_ploidy": 4, "agreement": "agree"},
}


@pytest.mark.parametrize(
    ("changes", "expected_reason"),
    [
        pytest.param({"level": None}, None, id="No level"),
        pytest.param({"status": "from_reference", "confidence": 0.1}, None, id="Not paper-derived"),
        pytest.param({}, None, id="Sound call"),
        pytest.param(
            {"formula_backed": True, "anchored": False, "confidence": 0.1}, None, id="Formula-backed"
        ),
        pytest.param({"method": "consensus_3way", "anchored": False}, None, id="3-way consensus"),
        pytest.param({"method": "rule_wins (conflict)"}, "conflicting_evidence", id="Methods conflict"),
        pytest.param({"anchored": False}, "not_anchored_to_target", id="Not anchored"),
        pytest.param({"level": 2, "haploid_mentions": 2}, "haploid_counter_evidence", id="Haploid mentions"),
        pytest.param({"haploid_mentions": 2}, None, id="Haploid mentions only matter for diploid"),
        pytest.param({"confidence": 0.49}, "weak_evidence", id="Low confidence"),
        pytest.param(
            {"reference_check": {"goat_ploidy": 2, "agreement": "disagree"}},
            "disagrees_with_reference",
            id="Disagrees with GoaT",
        ),
    ],
)
def test_precision_guard_reasons(changes: dict, expected_reason: str | None) -> None:
    """Test :func:`pipeline._precision_guard` keeps sound calls and withholds the rest, with the reason.

    Args:
        changes: Changes to a sound paper-derived ploidy result.
        expected_reason: Expected withholding reason, or None if the result is kept.

    """
    result = pipeline._precision_guard({**SOUND, **changes}, {"reference_ploidy": None})
    assert result.get("withheld_paper_ploidy", {}).get("reason") == expected_reason


@pytest.mark.parametrize(
    ("changes", "reference_ploidy", "expected"),
    [
        pytest.param({"anchored": False}, 6, (6, "from_reference", "GoaT"), id="Replaced by GoaT"),
        pytest.param({"anchored": False}, None, (None, "withheld", None), id="Withheld without GoaT"),
        pytest.param(
            {"priority_derivation": "doubled_monoploid", "formula_backed": True},
            2,
            (None, "withheld", None),
            id="Doubled monoploid is not diploid, GoaT not used",
        ),
    ],
)
def test_precision_guard_output(changes: dict, reference_ploidy: int | None, expected: tuple) -> None:
    """Test :func:`pipeline._precision_guard` replaces a withheld level with GoaT's when it may.

    Args:
        changes: Changes to a sound paper-derived ploidy result.
        reference_ploidy: GoaT reference ploidy of the assembly.
        expected: Expected ``(level, status, source)``.

    """
    result = pipeline._precision_guard({**SOUND, **changes}, {"reference_ploidy": reference_ploidy})
    assert (result["level"], result["status"], result["source"]) == expected
    assert result["withheld_paper_ploidy"]["level"] == 4


@pytest.mark.parametrize(
    ("result", "expected"),
    [
        pytest.param({"ploidy": {"level": 4}}, True, id="Level recovered"),
        pytest.param({"ploidy": {"level": None}}, False, id="No level"),
        pytest.param({"error": "boom", "ploidy": {"level": 4}}, False, id="Error"),
        pytest.param({}, False, id="No ploidy"),
    ],
)
def test_ploidy_recovered(result: dict, expected: bool) -> None:
    """Test :func:`pipeline._ploidy_recovered`.

    Args:
        result: Pipeline result.
        expected: Whether a ploidy level was recovered.

    """
    assert pipeline._ploidy_recovered(result) == expected


@pytest.mark.parametrize(
    ("paper", "accession", "expected"),
    [
        pytest.param({"retrieval_source": "linked_pmid"}, ACCESSION, True, id="Linked PMID"),
        pytest.param({"retrieval_source": "elink"}, ACCESSION, True, id="Entrez elink"),
        pytest.param({"retrieval_source": "bioproject"}, ACCESSION, True, id="BioProject"),
        pytest.param({"is_reference_paper": True}, ACCESSION, True, id="Reference paper"),
        pytest.param(
            {"text_data": {"text": "Data: gca_003086295.1"}}, ACCESSION, True, id="Accession in text"
        ),
        pytest.param(
            {"text_data": {"text": "Arachis hypogaea is tetraploid."}}, ACCESSION, False, id="Name only"
        ),
        pytest.param({"text_data": {"text": "anything"}}, "", False, id="No accession"),
    ],
)
def test_is_correct_paper_for_accession(paper: dict, accession: str, expected: bool) -> None:
    """Test :func:`pipeline._is_correct_paper_for_accession`.

    Args:
        paper: Paper record.
        accession: Assembly accession.
        expected: Whether Gemma may read the paper.

    """
    assert pipeline._is_correct_paper_for_accession(paper, {"assembly_accession": accession}) == expected


@pytest.mark.parametrize(
    ("assembly", "papers", "n_sections", "expected"),
    [
        pytest.param(
            {"assembly_accession": ACCESSION},
            [],
            1,
            ("error", "no assembly metadata"),
            id="No assembly metadata",
        ),
        pytest.param(ASSEMBLY, [], 1, ("success", None), id="No paper, GoaT used"),
        pytest.param(
            {**ASSEMBLY, "reference_ploidy": None},
            [],
            1,
            ("error", "No relevant papers found"),
            id="No paper, no GoaT",
        ),
        pytest.param(ASSEMBLY, [_paper()], 0, ("error", "no sections extracted"), id="Nothing parsed"),
    ],
)
@patch("ensembl.io.genomio.literature.pipeline.run_gemma_ploidy")
@patch("ensembl.io.genomio.literature.pipeline.parse_paper")
@patch("ensembl.io.genomio.literature.pipeline.fetch_papers_for_assembly")
def test_run_pipeline_early_exits(
    mock_fetch: MagicMock,
    mock_parse: MagicMock,
    mock_gemma: MagicMock,
    assembly: dict,
    papers: list[dict],
    n_sections: int,
    expected: tuple,
) -> None:
    """Test :func:`pipeline.run_pipeline` stops early, without Gemma, when there is nothing to read.

    Args:
        mock_fetch: Patched `fetch_papers_for_assembly()`.
        mock_parse: Patched `parse_paper()`.
        mock_gemma: Patched `run_gemma_ploidy()`.
        assembly: Assembly metadata.
        papers: Accepted papers.
        n_sections: Number of sections parsed from the best paper.
        expected: Expected ``(status, error)``.

    """
    mock_fetch.return_value = {"assembly": assembly, "papers": papers}
    mock_parse.return_value = {"n_sections": n_sections}
    result = pipeline.run_pipeline(ACCESSION)
    assert (result["status"], result.get("error")) == expected
    mock_gemma.assert_not_called()


@pytest.mark.parametrize(
    ("ensemble_ploidy", "gemma_on", "paper", "gemma_called", "expected"),
    [
        pytest.param(dict(SOUND), False, _paper(), False, (4, "success"), id="Paper ploidy kept"),
        pytest.param(
            {**SOUND, "level": 6, "anchored": False},
            False,
            _paper(),
            False,
            (4, "success"),
            id="Weak, GoaT used",
        ),
        pytest.param({**SOUND, "level": None}, False, _paper(), False, (None, "error"), id="No level"),
        pytest.param(dict(SOUND), True, _paper(), True, (4, "success"), id="Gemma reads the tied paper"),
        pytest.param(
            dict(SOUND),
            True,
            _paper("europepmc_name"),
            False,
            (4, "success"),
            id="Gemma skips an untied paper",
        ),
    ],
)
@patch(
    "ensembl.io.genomio.literature.pipeline.combine_with_gemma", lambda ensemble, gemma: {**ensemble, **gemma}
)
@patch("ensembl.io.genomio.literature.pipeline.extract_metadata", MagicMock())
@patch("ensembl.io.genomio.literature.pipeline.run_gemma_ploidy", MagicMock(return_value={"gemma": "read"}))
@patch("ensembl.io.genomio.literature.pipeline.is_enabled")
@patch("ensembl.io.genomio.literature.pipeline.run_vector_search")
@patch("ensembl.io.genomio.literature.pipeline.parse_paper")
@patch("ensembl.io.genomio.literature.pipeline.fetch_papers_for_assembly")
def test_run_pipeline(
    mock_fetch: MagicMock,
    mock_parse: MagicMock,
    mock_vector: MagicMock,
    mock_enabled: MagicMock,
    ensemble_ploidy: dict,
    gemma_on: bool,
    paper: dict,
    gemma_called: bool,
    expected: tuple,
) -> None:
    """Test :func:`pipeline.run_pipeline` on the main path: Gemma gate, GoaT check and precision guard.

    Args:
        mock_fetch: Patched `fetch_papers_for_assembly()`.
        mock_parse: Patched `parse_paper()`.
        mock_vector: Patched `run_vector_search()`.
        mock_enabled: Patched `is_enabled()`.
        ensemble_ploidy: Ploidy result of the rule + vector ensemble.
        gemma_on: Whether the Gemma layer is enabled.
        paper: The accepted paper.
        gemma_called: Whether Gemma is expected to read the paper.
        expected: Expected ``(ploidy level, status)``.

    """
    mock_fetch.return_value = {"assembly": ASSEMBLY, "papers": [paper], "reference_paper": {"title": "Ref"}}
    mock_parse.return_value = {"n_sections": 1, "n_chunks": 1, "sections": {"Abstract": "a"}}
    mock_vector.return_value = {"species": {"value": "Arachis hypogaea"}, "ploidy": ensemble_ploidy}
    mock_enabled.return_value = gemma_on
    final = pipeline.run_pipeline(ACCESSION)
    assert (final["ploidy"]["level"], final["status"]) == expected
    assert (final["paper"]["title"], final["reference_paper"]) == ("Peanut genome", {"title": "Ref"})
    assert (final["ploidy"].get("gemma") == "read") == gemma_called


@pytest.mark.xfail(strict=True, reason="pipeline.py:357 passes only final['ploidy'] to combine_with_gemma")
@patch("ensembl.io.genomio.literature.pipeline.is_enabled", lambda: True)
@patch("ensembl.io.genomio.literature.pipeline.extract_metadata", MagicMock())
@patch("ensembl.io.genomio.literature.pipeline.run_gemma_ploidy")
@patch("ensembl.io.genomio.literature.pipeline.run_vector_search")
@patch("ensembl.io.genomio.literature.pipeline.parse_paper")
@patch("ensembl.io.genomio.literature.pipeline.fetch_papers_for_assembly")
def test_run_pipeline_gemma_backfills_sex_and_cultivar(
    mock_fetch: MagicMock, mock_parse: MagicMock, mock_vector: MagicMock, mock_gemma: MagicMock
) -> None:
    """Test :func:`pipeline.run_pipeline` puts Gemma's sex and cultivar into the final record.

    `combine_with_gemma` is documented to backfill sex and cultivar when they are missing.

    Args:
        mock_fetch: Patched `fetch_papers_for_assembly()`.
        mock_parse: Patched `parse_paper()`.
        mock_vector: Patched `run_vector_search()`.
        mock_gemma: Patched `run_gemma_ploidy()`.

    """
    mock_fetch.return_value = {"assembly": ASSEMBLY, "papers": [_paper()]}
    mock_parse.return_value = {"n_sections": 1, "n_chunks": 1, "sections": {"Abstract": "a"}}
    mock_vector.return_value = {
        "species": {"value": "Arachis hypogaea"},
        "ploidy": dict(SOUND),
        "sex": {"value": "unknown"},
        "cultivars": {"confirmed": [], "unconfirmed": []},
    }
    mock_gemma.return_value = {"level": 4, "sex": "female", "strain_cultivar": "Tifrunner"}
    final = pipeline.run_pipeline(ACCESSION)
    assert (final["sex"]["value"], final["cultivars"]["unconfirmed"]) == ("female", ["Tifrunner"])


@patch("ensembl.io.genomio.literature.pipeline.run_pipeline")
def test_run_batch(mock_run: MagicMock) -> None:
    """Test :func:`pipeline.run_batch` keeps going after a failed accession.

    Args:
        mock_run: Patched `run_pipeline()`.

    """
    mock_run.side_effect = [{"assembly_accession": "A", "status": "success"}, RuntimeError("NCBI down")]
    assert pipeline.run_batch(["A", "B"], max_papers=2) == [
        {"assembly_accession": "A", "status": "success"},
        {"assembly_accession": "B", "status": "error", "error": "NCBI down"},
    ]
    assert mock_run.call_args_list == [call("A", max_papers=2), call("B", max_papers=2)]


@pytest.mark.parametrize(
    ("argv", "expectation"),
    [
        pytest.param([], pytest.raises(SystemExit), id="No source"),
        pytest.param(["--accession", "A", "--batch", "f"], pytest.raises(SystemExit), id="Both sources"),
        pytest.param(
            ["--accession", "A", "--max_papers", "x"], pytest.raises(SystemExit), id="Bad max_papers"
        ),
        pytest.param(["--accession", ACCESSION], does_not_raise(), id="Accession"),
    ],
)
@patch("ensembl.io.genomio.literature.pipeline.run_pipeline", return_value={"status": "success"})
def test_main_arguments(mock_run: MagicMock, argv: list[str], expectation: ContextManager) -> None:
    """Test :func:`pipeline.main` argument parsing.

    Args:
        mock_run: Patched `run_pipeline()`.
        argv: Command-line arguments.
        expectation: Context manager for the expected outcome.

    """
    with expectation:
        pipeline.main(argv)
        mock_run.assert_called_once_with(ACCESSION, max_papers=5)


@pytest.mark.parametrize(
    ("batch", "expected_call"),
    [
        pytest.param(False, call(ACCESSION, max_papers=2), id="Accession to stdout"),
        pytest.param(True, call(["GCA_1.1", "GCA_2.1"], max_papers=2), id="Batch to file, Gemma on"),
    ],
)
@patch.dict("os.environ")
@patch("ensembl.io.genomio.literature.pipeline.run_batch")
@patch("ensembl.io.genomio.literature.pipeline.run_pipeline")
def test_main(
    mock_run: MagicMock,
    mock_batch: MagicMock,
    batch: bool,
    expected_call: Any,
    tmp_path: Path,
    capsys: pytest.CaptureFixture,
) -> None:
    """Test :func:`pipeline.main` runs one accession or a batch and writes the JSON result.

    Args:
        mock_run: Patched `run_pipeline()`.
        mock_batch: Patched `run_batch()`.
        batch: Whether to run a batch file (written to a file, with ``--gemma``) or one accession.
        expected_call: Expected call to `run_batch()` or `run_pipeline()`.
        tmp_path: Pytest-provided temporary directory.
        capsys: Pytest fixture to capture stdout.

    """
    os.environ.pop("GEMMA_ENABLED", None)
    mock_run.return_value = mock_batch.return_value = {"species": "Arachis hypogaea × A. monticola"}
    batch_file = tmp_path / "accessions.txt"
    batch_file.write_text(" GCA_1.1 \n\nGCA_2.1\n")
    output = tmp_path / "out.json"
    if batch:
        pipeline.main(["--batch", str(batch_file), "--max_papers", "2", "--output", str(output), "--gemma"])
        written = output.read_text(encoding="utf-8")
    else:
        pipeline.main(["--accession", ACCESSION, "--max_papers", "2"])
        written = capsys.readouterr().out
    assert (mock_batch if batch else mock_run).call_args == expected_call
    assert json.loads(written) == mock_run.return_value
    assert "×" in written
    assert (os.environ.get("GEMMA_ENABLED") == "1") == batch


@pytest.mark.parametrize(
    ("name", "expected", "expectation"),
    [
        pytest.param("run_pipeline", pipeline.run_pipeline, does_not_raise(), id="run_pipeline"),
        pytest.param("run_batch", pipeline.run_batch, does_not_raise(), id="run_batch"),
        pytest.param(
            "no_such_name", None, pytest.raises(AttributeError, match="no_such_name"), id="Unknown name"
        ),
    ],
)
def test_package_lazy_attributes(name: str, expected: object, expectation: ContextManager) -> None:
    """Test the package exposes `run_pipeline` and `run_batch` lazily and nothing else.

    Args:
        name: Attribute requested from the package.
        expected: Expected attribute.
        expectation: Context manager for the expected exception, if any.

    """
    with expectation:
        assert getattr(literature, name) is expected
