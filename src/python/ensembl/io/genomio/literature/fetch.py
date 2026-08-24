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
"""Score candidate papers and orchestrate literature metadata retrieval for an assembly."""

import logging

import requests

from .metadata import fetch_assembly_metadata, enrich_assembly_metadata
from .ncbi import (
    fetch_linked_pubmed_for_assembly,
    fetch_bioproject_reference_papers,
    search_ncbi_entrez,
)
from .europepmc import (
    search_by_pmid,
    search_europe_pmc_by_name,
    search_europe_pmc_fulltext,
    get_fulltext_by_pmcid,
    validate_paper_ids,
    get_best_available_text,
)

logger = logging.getLogger(__name__)


ORGANELLE_KEYWORDS = [
    "chloroplast",
    "mitochondria",
    "plastid",
    "organelle",
    "cpdna",
    "mtdna",
    "plastome",
]


GENOME_KEYWORDS = [
    "genome assembly",
    "genome sequence",
    "reference genome",
    "chromosome-level",
    "chromosome-scale",
    "de novo assembly",
    "whole genome",
    "whole-genome",
    "genome of",
]


# Phrases that mark the target organism as a *habitat/host*, not the subject
# of the paper — e.g. "Draft genome of Bacillus ... isolated from cotton". Such
# papers mention the target's name only as the source and otherwise score high
# on name + genome keywords, so they need a hard penalty (see scorer).
HOST_CONTEXT_KEYWORDS = [
    "isolated from",
    "associated with",
    "endophyt",
    "rhizosphere",
    "phyllosphere",
    "symbiont",
    "symbiotic bacteri",
    "pathogen of",
    "microbiome",
    "microbiota",
    "gut of",
    "gut microb",
]


PLOIDY_HINTS = [
    "diploid",
    "triploid",
    "tetraploid",
    "pentaploid",
    "hexaploid",
    "octoploid",
    "polyploid",
    "allopolyploid",
    "autopolyploid",
    "amphidiploid",
    "ploidy",
    "2n=",
    "2n =",
]


# A candidate scoring >= STRONG_SCORE is treated as a confident genome paper,
# so Mode B will NOT escalate to extra name-based searches.
STRONG_SCORE = 5.0


def score_paper_candidate(paper: dict, identifiers: dict) -> float:
    """Relevance score for "is this the genome/ploidy paper for THIS species?".

    Used both to rank pooled candidates and to decide Mode-B escalation.
    """
    title = (paper.get("title") or "").lower()
    abstract = (paper.get("abstractText") or "").lower()
    text = f"{title} {abstract}"

    sci = (identifiers.get("scientific_name") or "").lower()
    common = (identifiers.get("common_name") or "").lower()
    genus = sci.split()[0] if sci else ""

    score = 0.0

    # name match (scientific > common > genus)
    if sci and sci in title:
        score += 3.0
    elif sci and sci in text:
        score += 1.5
    if common and len(common) >= 4 and common in title:
        score += 1.5
    if genus and genus in title:
        score += 1.0

    # genome-paper signal (word-order tolerant: "assembly of the ... genome")
    if any(k in title for k in GENOME_KEYWORDS):
        score += 3.0
    elif "genome" in title and "assembl" in title:
        score += 2.5
    elif "genome" in title or "assembl" in title:
        score += 1.5
    elif any(k in text for k in GENOME_KEYWORDS):
        score += 1.0

    # ploidy signal
    if any(k in text for k in PLOIDY_HINTS):
        score += 1.0

    # Full-text availability: a PMCID means the article's full text is open access
    # in Europe PMC, so the pipeline can read the whole body (where ploidy,
    # chromosome and cultivar statements usually sit) rather than only the
    # abstract. Papers with retrievable full text are therefore preferred.
    if paper.get("pmcid") or paper.get("pmCid"):
        score += 2.0

    # BioProject-linked reference paper: canonical publication for the assembly
    if paper.get("is_reference_paper"):
        score += 2.5

    # Full-text-confirmed species source: the exact binomial was found in the
    # article body (confirmation pass in fetch_papers_for_assembly) even though
    # it is absent from title/abstract — a multi-species / methods paper that
    # assembled this genome. Reward it like a strong name-in-title match so it
    # can outrank an off-target paper that merely names the species up front.
    if paper.get("_fulltext_name_confirmed"):
        score += 4.0

    # organelle penalty (a chloroplast/mito paper must never rank as a strong
    # nuclear-genome candidate, so the title penalty is large)
    if any(k in title for k in ORGANELLE_KEYWORDS):
        score -= 10.0
    elif any(k in abstract for k in ORGANELLE_KEYWORDS):
        score -= 2.0

    # host-context penalty: a paper about a microbe *isolated from* / *associated
    # with* the target (e.g. "Draft genome of Bacillus ... isolated from cotton")
    # names the target only as its habitat. These otherwise rank top on
    # name + genome-keyword hits, so penalise as hard as organelle papers.
    if any(k in title for k in HOST_CONTEXT_KEYWORDS):
        score -= 10.0
    elif any(k in abstract for k in HOST_CONTEXT_KEYWORDS):
        score -= 2.0

    return score


def _dedupe_and_rank(candidates: list, identifiers: dict) -> list:
    """Drop duplicate papers across sources, attach relevance_score, sort desc."""
    seen = set()
    rank = []
    for paper in candidates:
        key = (
            (paper.get("pmid") or "").strip()
            or (paper.get("pmcid") or paper.get("pmCid") or "").strip()
            or (paper.get("doi") or "").strip()
            or (paper.get("title", "")[:80].lower())
        )
        if not key or key in seen:
            continue
        seen.add(key)
        paper["relevance_score"] = round(score_paper_candidate(paper, identifiers), 2)
        rank.append(paper)
    rank.sort(key=lambda x: x["relevance_score"], reverse=True)
    return rank


def has_usable_fulltext(papers: list[dict]) -> bool:
    """Return True if at least one paper has a PMCID (open-access full text) and
    is not an organelle (chloroplast/mitochondrion) genome paper."""
    for p in papers:
        pmcid = p.get("pmcid") or p.get("pmCid") or ""
        if not pmcid:
            continue
        title = p.get("title", "").lower()
        abstract = (p.get("abstractText", "") or "").lower()
        is_organelle = any(kw in title or kw in abstract for kw in ORGANELLE_KEYWORDS)
        if not is_organelle:
            return True
    return False


def search_semantic_scholar(scientific_name: str, max_results: int = 5) -> list:
    """Last-resort search: query Semantic Scholar for genome-assembly papers by
    scientific name, keeping only genus-relevant hits. Returns Europe-PMC-shaped
    records ([] on failure)."""
    logger.info(f"  [Semantic Scholar] Searching for {scientific_name} ...")

    url = "https://api.semanticscholar.org/graph/v1/paper/search"
    try:
        response = requests.get(
            url,
            params={  # type: ignore[arg-type]
                "query": f"{scientific_name} genome assembly ploidy",
                "fields": "title,abstract,externalIds,openAccessPdf",
                "limit": max_results * 2,
            },
            timeout=30,
        )
        response.raise_for_status()
        papers_raw = response.json().get("data", [])
    except requests.RequestException as e:
        logger.warning(f"  [Semantic Scholar] Search failed: {e}")
        return []

    genus = scientific_name.split()[0].lower()
    papers = []
    for p in papers_raw:
        title = p.get("title", "")
        abstract = p.get("abstract", "") or ""
        if genus not in (title + abstract).lower():
            continue
        ext_ids = p.get("externalIds", {})
        papers.append(
            {
                "title": title,
                "pmcid": f"PMC{ext_ids['PubMedCentral']}" if "PubMedCentral" in ext_ids else "",
                "pmid": str(ext_ids.get("PubMed", "")),
                "doi": ext_ids.get("DOI", ""),
                "abstractText": abstract,
                "source": "semantic_scholar",
            }
        )
        logger.info(f"  [Semantic Scholar] {title[:60]}...")
        if len(papers) >= max_results:
            break

    logger.info(f"  [Semantic Scholar] Found {len(papers)} relevant paper(s)")
    return papers


def _strongest_score(candidates: list, identifiers: dict) -> float:
    """Highest relevance score among candidates (0.0 if the list is empty)."""
    return max((score_paper_candidate(p, identifiers) for p in candidates), default=0.0)


def _gather_paper_candidates(accession: str, assembly: dict, identifiers: dict, max_results: int) -> list:
    """Run the Mode-B retrieval escalation and return the pooled candidate papers.

    Sources are tried in priority order (directly linked PMIDs -> Entrez elink ->
    BioProject reference -> multi-identifier name search -> full-text species
    search -> Entrez PMC -> Semantic Scholar) and only escalated to the weaker
    sources while no STRONG genome paper has surfaced. Each paper is tagged with
    its `retrieval_source`."""
    scientific_name = identifiers["scientific_name"]
    common_name = identifiers["common_name"]
    taxon_id = identifiers["taxon_id"]
    candidates: list = []

    # Priority 1: directly linked PMIDs from the NCBI assembly record
    linked_pmids = assembly.get("linked_pmids", [])
    if linked_pmids:
        logger.info(f"\n[2/4] Fetching {len(linked_pmids)} directly linked paper(s) ...")
        for pmid in linked_pmids:
            paper = search_by_pmid(pmid)
            if paper:
                paper["retrieval_source"] = "linked_pmid"
                candidates.append(paper)
                logger.info(f"      Found: {paper.get('title', '')[:70]}...")

    # Priority 2: NCBI Entrez elink (curated assembly -> pubmed)
    if not candidates:
        logger.info(f"\n[2/4] Trying NCBI Entrez elink for {accession} ...")
        for p in fetch_linked_pubmed_for_assembly(accession, max_results=max_results):
            p["retrieval_source"] = "elink"
            candidates.append(p)

    # Priority 2b: BioProject reference paper (assembly -> bioproject -> pubmed).
    # Always pooled — the canonical reference paper should be attachable even
    # when it carries no ploidy information (mentor item 3).
    logger.info(f"\n[2/4] Looking up BioProject reference paper for {accession} ...")
    for p in fetch_bioproject_reference_papers(accession, max_results=max_results):
        candidates.append(p)

    # ---- Mode B escalation: if the accession-linked candidates contain no
    #      STRONG genome paper, search by scientific + common name and merge.
    if _strongest_score(candidates, identifiers) < STRONG_SCORE:
        logger.info(
            f"\n[2/4] No strong genome paper from accession links "
            f"(best score {_strongest_score(candidates, identifiers):.1f}) "
            f"-> multi-identifier name search ..."
        )
        for p in search_europe_pmc_by_name(taxon_id, scientific_name, common_name, max_results=max_results):
            p.setdefault("retrieval_source", "europepmc_name")
            candidates.append(p)

    # Priority 3c: full-text species search. When no clearly-strong genome paper
    # has surfaced, the true source may be a multi-species / methods paper that
    # names this species only in its body (e.g. an assembly reported alongside
    # others). Search Europe PMC full text and pool the open-access hits; a
    # confirmation pass below verifies each before it can rank.
    if _strongest_score(candidates, identifiers) < STRONG_SCORE + 2.0:
        for p in search_europe_pmc_fulltext(scientific_name, max_results=max_results):
            p.setdefault("retrieval_source", "europepmc_fulltext")
            candidates.append(p)

    # Priority 4: NCBI Entrez PMC search — still no strong candidate
    if _strongest_score(candidates, identifiers) < STRONG_SCORE:
        logger.info(f"\n[2/4] Trying NCBI Entrez PMC search for {scientific_name} ...")
        for p in search_ncbi_entrez(scientific_name, max_results=max_results):
            p.setdefault("retrieval_source", "entrez")
            candidates.append(p)

    # Priority 5: Semantic Scholar — nothing found at all
    if not candidates:
        logger.info(f"\n[2/4] Trying Semantic Scholar for {scientific_name} ...")
        for p in search_semantic_scholar(scientific_name, max_results=max_results):
            p.setdefault("retrieval_source", "semantic_scholar")
            candidates.append(p)

    if not candidates:
        logger.info(f"\n[2/4] All search sources exhausted. No papers found.")

    return candidates


def _confirm_fulltext_candidates(candidates: list, scientific_name: str) -> list:
    """Verify full-text-only species matches and drop spurious ones.

    Candidates tagged `_needs_fulltext_confirm` matched the species only in their
    full text. Fetch the text and keep only those where the exact binomial really
    appears (marking them so scoring can reward a genuine source); drop spurious
    matches. Bounded to a few fetches so cost stays predictable."""
    sci_lower = scientific_name.lower()
    confirmed, ft_checks = [], 0
    for p in candidates:
        if not p.get("_needs_fulltext_confirm"):
            confirmed.append(p)
            continue
        pmcid = p.get("pmcid") or p.get("pmCid")
        if pmcid and ft_checks < 5:
            ft_checks += 1
            ft = get_fulltext_by_pmcid(pmcid)
            if ft and sci_lower in ft.lower():
                p["_fulltext_name_confirmed"] = True
                p.pop("_needs_fulltext_confirm", None)
                confirmed.append(p)
                logger.info(
                    f"      [fulltext-confirm] {scientific_name} found in {pmcid}: "
                    f"{p.get('title', '')[:55]}..."
                )
            # else: species not actually in the body -> drop spurious match
    return confirmed


def _prepare_papers(papers_found: list) -> list:
    """Validate IDs and attach the best available text for each ranked paper."""
    logger.info(f"\n[3/4] Validating IDs and fetching text ...")
    enriched_papers = []
    for paper in papers_found:
        validation = validate_paper_ids(paper)
        text_data = get_best_available_text(paper)
        enriched_papers.append(
            {
                "title": paper.get("title", ""),
                "pmcid": paper.get("pmcid"),
                "pmid": paper.get("pmid"),
                "doi": paper.get("doi"),
                "retrieval_source": paper.get("retrieval_source"),
                "relevance_score": paper.get("relevance_score"),
                "is_reference_paper": bool(paper.get("is_reference_paper")),
                "validation": validation,
                "text_data": text_data,
            }
        )
        logger.info(
            f"      [{validation['best_available']:>10}]  "
            f"(score {paper.get('relevance_score')}, {paper.get('retrieval_source')})  "
            f"{paper.get('title', '')[:52]}..."
        )
    return enriched_papers


def _extract_reference_paper(papers_found: list) -> dict | None:
    """Return the canonical BioProject-linked reference paper, if one is present.

    Attached regardless of whether it ends up carrying ploidy (mentor item 3)."""
    for paper in papers_found:
        if paper.get("is_reference_paper"):
            return {
                "title": paper.get("title", ""),
                "pmcid": paper.get("pmcid"),
                "pmid": paper.get("pmid"),
                "doi": paper.get("doi"),
                "retrieval_source": "bioproject",
                "relevance_score": paper.get("relevance_score"),
            }
    return None


def fetch_papers_for_assembly(accession: str, max_results: int = 5) -> dict:
    """Retrieve assembly metadata and the ranked candidate papers for one accession.

    Orchestrates the workflow stages and returns {assembly, papers, reference_paper}:
      1. fetch + enrich assembly metadata (NCBI Datasets, Taxonomy, GoaT);
      2. gather candidate papers (`_gather_paper_candidates`);
      3. confirm full-text-only species matches (`_confirm_fulltext_candidates`);
      4. de-duplicate and rank (`_dedupe_and_rank`);
      5. validate IDs and fetch text (`_prepare_papers`)."""
    logger.info(f"\n[1/4] Fetching assembly metadata for {accession} ...")
    assembly = fetch_assembly_metadata(accession)
    assembly = enrich_assembly_metadata(assembly)  # common_name + GoaT reference via taxon_id
    # fetch_assembly_metadata may have resolved a stale version to the latest;
    # use that for all downstream accession-based lookups (elink, BioProject).
    accession = assembly.get("assembly_accession", accession) or accession
    scientific_name = assembly.get("scientific_name", "") or ""
    common_name = assembly.get("common_name", "") or ""
    taxon_id = assembly.get("taxon_id", "") or ""

    # Mentor's "key metadata fields" — used as multiple retrieval identifiers
    identifiers = {
        "accession": accession,
        "assembly_name": assembly.get("assembly_name"),
        "scientific_name": scientific_name,
        "common_name": common_name,
        "taxon_id": taxon_id,
    }

    logger.info(f"      taxon_id          = {taxon_id}")
    logger.info(f"      scientific_name   = {scientific_name}")
    logger.info(f"      common_name       = {common_name or '(none)'}")
    logger.info(f"      assembly_name     = {assembly.get('assembly_name')}")
    logger.info(
        f"      chromosome_number = {assembly.get('chromosome_number')} "
        f"(source: {assembly.get('chromosome_source') or 'not found -> will try text extraction'})"
    )

    # Guard: no scientific name / taxon -> name-based search would degenerate
    # into an empty query and return unrelated papers. Skip retrieval entirely.
    if not scientific_name and not taxon_id:
        logger.info("      No NCBI scientific_name / taxon_id — skipping paper retrieval.")
        return {"assembly": assembly, "papers": []}

    candidates = _gather_paper_candidates(accession, assembly, identifiers, max_results)
    candidates = _confirm_fulltext_candidates(candidates, scientific_name)

    # ---- de-duplicate across sources + rank best-first by relevance ----
    papers_found = _dedupe_and_rank(candidates, identifiers)
    if papers_found:
        top = papers_found[0]
        logger.info(
            f"\n      {len(papers_found)} unique candidate(s) ranked; "
            f"top score {top['relevance_score']} ({top.get('retrieval_source')})"
        )

    enriched_papers = _prepare_papers(papers_found)
    reference_paper = _extract_reference_paper(papers_found)

    logger.info(f"\n[4/4] Done. {len(enriched_papers)} paper(s) ready for extraction.")
    logger.info(
        f"      chromosome_number will be "
        f"{'used from NCBI' if assembly.get('chromosome_number') else 'extracted from text in extract.py'}"
    )

    return {
        "assembly": assembly,
        "papers": enriched_papers,
        "reference_paper": reference_paper,
    }
