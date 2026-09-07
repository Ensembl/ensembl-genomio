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
"""Orchestrate the literature metadata extraction pipeline and provide the CLI entry point."""

import argparse
import json
import logging

from ensembl.io.genomio.literature.extract import extract_metadata
from ensembl.io.genomio.literature.fetch import fetch_papers_for_assembly
from ensembl.io.genomio.literature.gemma import combine_with_gemma, is_enabled, run_gemma_ploidy
from ensembl.io.genomio.literature.parse import parse_paper
from ensembl.io.genomio.literature.search import run_vector_search

logger = logging.getLogger(__name__)

# ============================================================
# Paper quality filter
# ============================================================

EXCLUDE_KEYWORDS = [
    "chloroplast",
    "cpdna",
    "mitochondria",
    "mtdna",
    "organelle",
    "plastid",
    "plastome",
]


def is_nuclear_genome_paper(paper: dict) -> bool:
    """Return False if the paper is about an organelle (chloroplast / mitochondrion)
    genome rather than the nuclear genome, based on title and abstract keywords."""
    title = paper.get("title", "").lower()
    text = paper.get("text_data", {}).get("text") or ""
    abstract = text[:500].lower()  # check first 500 chars only

    return not any(kw in title or kw in abstract for kw in EXCLUDE_KEYWORDS)


def select_best_paper(papers: list[dict]) -> dict | None:
    """Pick the highest-relevance paper that has usable text, preferring full text
    over abstract on ties; falls back to the top candidate, or None if empty.

    Candidates arrive pre-sorted by relevance_score (desc) from fetch.py, whose
    scorer already rewards genome/ploidy/name signals and heavily penalises
    organelle papers, so relevance_score is a better selector than
    "any fulltext first".
    """
    usable = [
        paper
        for paper in papers
        if paper.get("text_data", {}).get("source") in ("fulltext", "abstract", "search_result")
    ]
    if not usable:
        if papers:
            logger.info(f"  Selected (fallback, no usable text): {papers[0].get('title', '')[:70]}...")
            return papers[0]
        return None

    # Rank by relevance_score; tie-break prefers fulltext over abstract.
    best = max(
        usable,
        key=lambda paper: (
            paper.get("relevance_score", 0.0),
            1 if paper.get("text_data", {}).get("source") == "fulltext" else 0,
        ),
    )
    logger.info(
        f"  Selected (score {best.get('relevance_score')}, "
        f"{best.get('retrieval_source')}, {best['text_data']['source']}): "
        f"{best.get('title', '')[:60]}..."
    )
    return best


def _goat_agreement(paper_level: int | None, goat_ploidy: int | None) -> dict:
    """Compare the paper-derived ploidy level against GoaT's species-level reference
    ploidy and return a curator-facing agreement flag (agree / disagree /
    no_reference / no_paper_ploidy). GoaT is species-level, so this flags conflicts
    rather than overriding the paper-derived value."""
    if goat_ploidy is None:
        agreement = "no_reference"
    elif paper_level is None:
        agreement = "no_paper_ploidy"
    elif paper_level == goat_ploidy:
        agreement = "agree"
    else:
        agreement = "disagree"
    return {"goat_ploidy": goat_ploidy, "agreement": agreement}


def _ploidy_recovered(result: dict) -> bool:
    """A batch result counts as useful only if a ploidy level was actually
    recovered: a paper that parsed fine but states no ploidy is not a success."""
    if result.get("error"):
        return False
    return (result.get("ploidy") or {}).get("level") is not None


def _is_correct_paper_for_accession(paper: dict, assembly: dict) -> bool:
    """Accession-level check before letting Gemma read a paper: only trust a paper
    that is genuinely tied to THIS assembly — linked by accession (directly linked
    PMID / Entrez elink) or by its BioProject, or with the assembly's scientific
    name actually present in the paper text. This stops Gemma from confidently
    reading ploidy out of an off-target paper."""
    if paper.get("retrieval_source") in ("linked_pmid", "elink", "bioproject") or paper.get(
        "is_reference_paper"
    ):
        return True
    scientific_name = (assembly.get("scientific_name") or "").lower()
    paper_text = (paper.get("text_data", {}).get("text") or "").lower()
    return bool(scientific_name) and scientific_name in paper_text


# ============================================================
# MAIN FUNCTION — run full pipeline for one assembly accession
#           input:  assembly accession string
#           output: final metadata dict
# ============================================================


def run_pipeline(accession: str, max_papers: int = 5) -> dict:
    """Run the full extraction pipeline for one assembly accession: fetch papers,
    select the best one, parse it, run rule-based extraction, then vector search
    and ensemble (plus the optional Gemma layer). Returns the final metadata dict."""
    logger.info("=" * 60)
    logger.info(f"Pipeline started for: {accession}")
    logger.info("=" * 60)

    # ── Step 1: fetch ─────────────────────────────────────────
    logger.info("\n[Step 1/4] Fetching assembly metadata and papers ...")
    fetched = fetch_papers_for_assembly(accession, max_results=max_papers)
    assembly = fetched["assembly"]
    papers = fetched["papers"]
    candidate_papers = fetched.get("candidate_papers", [])
    reference_paper = fetched.get("reference_paper")

    # If the strict gate accepted no accession/BioProject-linked paper, fall back to
    # the best candidate paper (found by name / full text) and treat its metadata as a
    # SUGGESTION rather than a confirmed result. This keeps a literature-grounded answer
    # (Gemma reads the whole paper; rule + vector still run) instead of dropping straight
    # to GoaT. Confirmed only when a paper was strictly linked to the accession.
    suggested = False
    if not papers and candidate_papers:
        logger.info("  No accession-linked paper; using the best candidate as a suggestion.")
        papers = candidate_papers
        suggested = True

    # Guard: NCBI returned no usable assembly metadata (stale / suppressed
    # accession). Without a scientific name we cannot search by name or trust
    # text-based species extraction, so stop here with a clear failure mode
    # instead of degenerating into an empty-query search and garbage output.
    if not assembly.get("scientific_name") and not assembly.get("taxon_id"):
        logger.info(
            "  No NCBI assembly metadata (scientific_name / taxon_id missing). "
            "Likely a stale or suppressed accession — stopping."
        )
        return {
            "assembly_accession": assembly.get("assembly_accession"),
            "status": "error",
            "error": "no assembly metadata",
        }

    if not papers:
        # Strict retrieval found no accession/BioProject-linked paper. Fall back to
        # GoaT's species-level reference ploidy when available, so we still return a
        # result (mentor: prioritise precision in literature retrieval, but use GoaT
        # where present); otherwise report nothing found.
        goat_level = assembly.get("reference_ploidy")
        if goat_level is not None:
            logger.info(f"  No qualifying paper — using GoaT reference ploidy ({goat_level}).")
            ploidy = {
                "level": goat_level,
                "mechanism": None,
                "status": "from_reference",
                "source": "GoaT",
                "method": "reference_taxon",
                "confidence": 0.6,
                "reference_check": {"goat_ploidy": goat_level, "agreement": "no_qualifying_paper"},
            }
            status, error = "success", None
        else:
            logger.info("  No qualifying paper and no GoaT reference — nothing to report.")
            ploidy = {"level": None, "status": "no_relevant_paper", "source": None}
            status, error = "error", "no relevant paper"

        result = {
            "assembly_accession": assembly.get("assembly_accession"),
            "assembly_name": assembly.get("assembly_name"),
            "taxon_id": assembly.get("taxon_id"),
            "scientific_name": assembly.get("scientific_name"),
            "chromosome_number": {
                "value": assembly.get("chromosome_number"),
                "source": assembly.get("chromosome_source"),
            },
            "ploidy": ploidy,
            "cultivars": {"confirmed": [], "unconfirmed": []},
            "sex": {"value": "unknown"},
            "paper": {"note": "No accession/BioProject-linked paper found"},
            "papers_found": 0,
            "status": status,
        }
        if error:
            result["error"] = error
        return result

    # Select best paper — prefer nuclear genome + fulltext
    logger.info(f"\n  Selecting best paper from {len(papers)} candidate(s) ...")
    best_paper = select_best_paper(papers)

    if best_paper is None:
        logger.info("  No usable paper found.")
        return {
            "assembly_accession": assembly.get("assembly_accession"),
            "status": "error",
            "error": "no usable paper",
            "paper": {"note": "Not found relevant paper"},
        }

    logger.info(f"  Text source: {best_paper['text_data']['source']}")
    logger.info(
        f"  Validated:   "
        f"PMCID={best_paper['validation']['pmcid_valid']} "
        f"PMID={best_paper['validation']['pmid_valid']}"
    )

    # ── Step 2: parse ─────────────────────────────────────────
    logger.info("\n[Step 2/4] Parsing paper text into sections ...")
    parsed = parse_paper(best_paper["text_data"])

    if parsed["n_sections"] == 0:
        logger.info("  No sections extracted. Check text source.")
        return {
            "assembly_accession": assembly.get("assembly_accession"),
            "assembly_name": assembly.get("assembly_name"),
            "taxon_id": assembly.get("taxon_id"),
            "status": "error",
            "error": "no sections extracted",
        }

    logger.info(f"  Sections: {parsed['n_sections']}  Chunks: {parsed['n_chunks']}")

    # ── Step 3: rule-based extraction ─────────────────────────
    logger.info("\n[Step 3/4] Running rule-based extraction ...")
    rule_result = extract_metadata(parsed, assembly, paper_title=best_paper.get("title"))

    # ── Step 4: vector search + ensemble ──────────────────────
    logger.info("\n[Step 4/4] Running vector search and ensemble ...")
    final = run_vector_search(parsed, rule_result)

    # ── Step 4.5 (optional): Gemma hybrid layer ───────────────
    #     Local Gemma reads the FULL paper (all parsed sections), not just the
    #     retrieved vector chunks, so it forms an independent opinion even when the
    #     decisive ploidy sentence didn't make the vector-search cut. combine_with_gemma
    #     folds its result in as a third vote. No-op unless GEMMA_ENABLED=1.
    if is_enabled() and not _is_correct_paper_for_accession(best_paper, assembly):
        logger.info(
            "  [gemma] skipped: the selected paper is not confidently tied to this "
            "accession (no accession/BioProject link and scientific name absent from text)"
        )
    elif is_enabled():
        gemma_passages = [{"text": text, "section": section} for section, text in parsed["sections"].items()]
        gemma_result = run_gemma_ploidy(
            final["species"]["value"],
            gemma_passages,
            accession=accession,
        )
        if gemma_result is not None:
            final["ploidy"] = combine_with_gemma(final["ploidy"], gemma_result)
            logger.info(
                f"  [gemma] combined -> level={final['ploidy']['level']} "
                f"(method: {final['ploidy']['method']})"
            )

    # ── GoaT cross-check ──────────────────────────────────────
    #     Flag (do not override) when the paper-derived ploidy disagrees with
    #     GoaT's species-level reference, so a curator can spot the conflict.
    final["ploidy"]["reference_check"] = _goat_agreement(
        final["ploidy"].get("level"), assembly.get("reference_ploidy")
    )
    if final["ploidy"]["reference_check"]["agreement"] == "disagree":
        logger.info(
            f"  [validate] paper ploidy {final['ploidy'].get('level')} disagrees with "
            f"GoaT reference {assembly.get('reference_ploidy')}"
        )

    # Add paper-level metadata to final output
    final["paper"] = {
        "title": best_paper["title"],
        "pmcid": best_paper["pmcid"],
        "pmid": best_paper["pmid"],
        "doi": best_paper["doi"],
        "text_source": best_paper["text_data"]["source"],
        "retrieval_source": best_paper.get("retrieval_source"),
        "relevance_score": best_paper.get("relevance_score"),
        "validation": best_paper["validation"],
    }
    final["papers_found"] = len(papers)
    final["reference_paper"] = reference_paper  # BioProject canonical paper (attached even if no ploidy)

    # When no paper was strictly linked to the accession, the answer came from the
    # best candidate paper — mark it as a suggestion (not a confirmed result) and
    # temper the confidence, so curators can tell literature-confirmed from inferred.
    if suggested and final["ploidy"].get("level") is not None:
        final["ploidy"]["status"] = "suggested"
        final["ploidy"]["source"] = "suggested_" + str(final["ploidy"].get("source") or "literature")
        final["ploidy"]["confidence"] = round(final["ploidy"].get("confidence", 0.5) * 0.7, 2)

    # A result is only "success" if a ploidy level was actually recovered; a paper
    # that parsed fine but states no ploidy is flagged as an error, so single-accession
    # (Nextflow) runs surface "nothing useful found" cleanly rather than as success.
    final["status"] = "success" if _ploidy_recovered(final) else "error"
    if final["status"] == "error":
        final.setdefault("error", "no ploidy recovered")

    return final


# ============================================================
# BATCH MODE — run pipeline for multiple accessions
#           input:  list of accession strings
#           output: list of final metadata dicts
# ============================================================


def run_batch(accessions: list[str], max_papers: int = 5) -> list[dict]:
    """Run run_pipeline over a list of accessions, tagging each result with a
    success/error status and continuing past failures. Returns the list of results."""
    results = []
    total = len(accessions)

    for i, accession in enumerate(accessions, 1):
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Batch progress: {i}/{total}")
        try:
            # run_pipeline already stamps result["status"] based on whether a
            # ploidy level was recovered; batch mode just handles hard failures.
            result = run_pipeline(accession, max_papers=max_papers)
        except Exception as e:
            logger.info(f"  ERROR for {accession}: {e}")
            result = {
                "assembly_accession": accession,
                "status": "error",
                "error": str(e),
            }
        results.append(result)

    logger.info(f"\n{'=' * 60}")
    logger.info(f"Batch complete: {len(results)} accessions processed")
    successful = sum(r.get("status") == "success" for r in results)
    logger.info(f"  Successful: {successful}/{total}")
    logger.info(f"  Failed:     {total - successful}/{total}")

    return results


# ============================================================
# CLI entry point
# ============================================================


def main(argv: list[str] | None = None) -> None:
    """CLI entry point: parse arguments and run a single accession (--accession) or
    a batch (--batch), optionally writing JSON results to --output (else stdout)."""
    # Surface fetch.py's logging output on the console (bare format matches the
    # previous print-based progress); the library itself stays log-config-free.
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    parser = argparse.ArgumentParser(description="Extract genomic metadata from an assembly accession.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--accession",
        help="Single NCBI assembly accession (e.g. GCA_003086295.3)",
    )
    source.add_argument(
        "--batch",
        help="Path to a text file with one accession per line (batch mode)",
    )
    parser.add_argument(
        "--max_papers",
        type=int,
        default=5,
        help="Max number of papers to fetch per accession (default: 5)",
    )
    parser.add_argument(
        "--output",
        help="Path to save results as JSON; printed to stdout if omitted",
    )
    args = parser.parse_args(argv)

    if args.accession:
        result = run_pipeline(args.accession, max_papers=args.max_papers)
        if args.output:
            with open(args.output, "w") as handle:
                json.dump(result, handle, indent=2, ensure_ascii=False)
            logger.info(f"\nResults saved to {args.output}")
        else:
            print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        with open(args.batch) as handle:
            accessions = [line.strip() for line in handle if line.strip()]
        logger.info(f"Batch mode: {len(accessions)} accessions loaded from {args.batch}")
        results = run_batch(accessions, max_papers=args.max_papers)
        if args.output:
            with open(args.output, "w") as handle:
                json.dump(results, handle, indent=2, ensure_ascii=False)
            logger.info(f"\nResults saved to {args.output}")
        else:
            print(json.dumps(results, indent=2, ensure_ascii=False))
