import argparse
import json

from .fetch import fetch_papers_for_assembly
from .parse import parse_paper
from .extract import extract_metadata
from .search import run_vector_search
from . import gemma
from .gemma import combine_with_gemma

# ============================================================
# Paper quality filter
# ============================================================

EXCLUDE_KEYWORDS = [
    "chloroplast", "mitochondria", "organelle",
    "plastid", "cpdna", "mtdna", "plastome",
]

def is_nuclear_genome_paper(paper: dict) -> bool:
    title    = paper.get("title", "").lower()
    text     = paper.get("text_data", {}).get("text") or ""
    abstract = text[:500].lower()  # check first 500 chars only

    return not any(
        kw in title or kw in abstract
        for kw in EXCLUDE_KEYWORDS
    )


def select_best_paper(papers: list[dict]) -> dict | None:
    # Candidates arrive pre-sorted by relevance_score (desc) from fetch.py.
    # We pick the highest-relevance paper that has usable text. relevance_score
    # already rewards genome/ploidy/name signals and heavily penalises organelle
    # papers, so it is a better selector than "any fulltext first".
    usable = [
        p for p in papers
        if p.get("text_data", {}).get("source") in ("fulltext", "abstract", "search_result")
    ]
    if not usable:
        if papers:
            print(f"  Selected (fallback, no usable text): {papers[0].get('title', '')[:70]}...")
            return papers[0]
        return None

    # Rank by relevance_score; tie-break prefers fulltext over abstract.
    def rank(p):
        return (
            p.get("relevance_score", 0.0),
            1 if p.get("text_data", {}).get("source") == "fulltext" else 0,
        )

    best = max(usable, key=rank)
    print(f"  Selected (score {best.get('relevance_score')}, "
          f"{best.get('retrieval_source')}, {best['text_data']['source']}): "
          f"{best.get('title', '')[:60]}...")
    return best


# ============================================================
# MAIN FUNCTION — run full pipeline for one assembly accession
#           input:  assembly accession string
#           output: final metadata dict
# ============================================================

def run_pipeline(accession: str, max_results: int = 5) -> dict:
    print("=" * 60)
    print(f"Pipeline started for: {accession}")
    print("=" * 60)

    # ── Step 1: fetch ─────────────────────────────────────────
    print("\n[Step 1/4] Fetching assembly metadata and papers ...")
    fetched  = fetch_papers_for_assembly(accession, max_results=max_results)
    assembly = fetched["assembly"]
    papers   = fetched["papers"]
    reference_paper = fetched.get("reference_paper")

    # Guard: NCBI returned no usable assembly metadata (stale / suppressed
    # accession). Without a scientific name we cannot search by name or trust
    # text-based species extraction, so stop here with a clear failure mode
    # instead of degenerating into an empty-query search and garbage output.
    if not assembly.get("scientific_name") and not assembly.get("taxon_id"):
        print("  No NCBI assembly metadata (scientific_name / taxon_id missing). "
              "Likely a stale or suppressed accession — stopping.")
        return {
            "assembly_accession": assembly.get("assembly_accession"),
            "error":              "no_assembly_metadata",
        }

    if not papers:
        print("  No relevant paper found. Returning assembly metadata only.")
        return {
            "assembly_accession": assembly.get("assembly_accession"),
            "assembly_name":      assembly.get("assembly_name"),
            "taxon_id":           assembly.get("taxon_id"),
            "scientific_name":    assembly.get("scientific_name"),
            "chromosome_number":  {
                "value":  assembly.get("chromosome_number"),
                "source": assembly.get("chromosome_source"),
            },
            "ploidy":    {"level": None, "status": "no_relevant_paper", "source": None},
            "cultivars": {"confirmed": [], "unconfirmed": []},
            "sex":       {"value": "unknown"},
            "paper":     {"note": "Not found relevant paper"},
            "papers_found": 0,
        }

    # Select best paper — prefer nuclear genome + fulltext
    print(f"\n  Selecting best paper from {len(papers)} candidate(s) ...")
    best_paper = select_best_paper(papers)

    if best_paper is None:
        print("  No usable paper found.")
        return {
            "assembly_accession": assembly.get("assembly_accession"),
            "error":              "no_usable_paper",
            "paper":              {"note": "Not found relevant paper"},
        }

    print(f"  Text source: {best_paper['text_data']['source']}")
    print(f"  Validated:   "
          f"PMCID={best_paper['validation']['pmcid_valid']} "
          f"PMID={best_paper['validation']['pmid_valid']}")

    # ── Step 2: parse ─────────────────────────────────────────
    print("\n[Step 2/4] Parsing paper text into sections ...")
    parsed = parse_paper(best_paper["text_data"])

    if parsed["n_sections"] == 0:
        print("  No sections extracted. Check text source.")
        return {
            "assembly_accession": assembly.get("assembly_accession"),
            "assembly_name":      assembly.get("assembly_name"),
            "taxon_id":           assembly.get("taxon_id"),
            "error":              "no_sections_extracted",
        }

    print(f"  Sections: {parsed['n_sections']}  Chunks: {parsed['n_chunks']}")

    # ── Step 3: rule-based extraction ─────────────────────────
    print("\n[Step 3/4] Running rule-based extraction ...")
    rule_result = extract_metadata(parsed, assembly, paper_title=best_paper.get("title"))

    # ── Step 4: vector search + ensemble ──────────────────────
    print("\n[Step 4/4] Running vector search and ensemble ...")
    final = run_vector_search(parsed, rule_result)

    # ── Step 4.5 (optional): Gemma hybrid layer ───────────────
    #     Local Gemma extracts a grounded ploidy candidate; combine_with_gemma
    #     folds it in as a third vote. No-op unless GEMMA_ENABLED=1.
    if gemma.is_enabled():
        g = gemma.run_gemma_ploidy(
            final["species"]["value"],
            final["ploidy"].get("evidence_passages", []),
            accession=accession,
        )
        if g is not None:
            final["ploidy"] = combine_with_gemma(final["ploidy"], g)
            print(f"  [gemma] combined -> level={final['ploidy']['level']} "
                  f"(method: {final['ploidy']['method']})")

    # Add paper-level metadata to final output
    final["paper"] = {
        "title":            best_paper["title"],
        "pmcid":            best_paper["pmcid"],
        "pmid":             best_paper["pmid"],
        "doi":              best_paper["doi"],
        "text_source":      best_paper["text_data"]["source"],
        "retrieval_source": best_paper.get("retrieval_source"),
        "relevance_score":  best_paper.get("relevance_score"),
        "validation":       best_paper["validation"],
    }
    final["papers_found"] = len(papers)
    final["reference_paper"] = reference_paper   # BioProject canonical paper (attached even if no ploidy)

    return final


# ============================================================
# BATCH MODE — run pipeline for multiple accessions
#           input:  list of accession strings
#           output: list of final metadata dicts
# ============================================================

def run_batch(accessions: list[str], max_results: int = 5) -> list[dict]:
    results = []
    total   = len(accessions)

    for i, accession in enumerate(accessions, 1):
        print(f"\n{'=' * 60}")
        print(f"Batch progress: {i}/{total}")
        try:
            result         = run_pipeline(accession, max_results=max_results)
            result["status"] = "success"
        except Exception as e:
            print(f"  ERROR for {accession}: {e}")
            result = {
                "assembly_accession": accession,
                "status":             "error",
                "error":              str(e),
            }
        results.append(result)

    print(f"\n{'=' * 60}")
    print(f"Batch complete: {len(results)} accessions processed")
    successful = sum(1 for r in results if r.get("status") == "success")
    print(f"  Successful: {successful}/{total}")
    print(f"  Failed:     {total - successful}/{total}")

    return results


# ============================================================
# CLI entry point
# ============================================================

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Extract genomic metadata from an assembly accession."
    )
    parser.add_argument(
        "--accession",
        type=str,
        help="Single NCBI assembly accession (e.g. GCA_003086295.3)",
    )
    parser.add_argument(
        "--batch",
        type=str,
        help="Path to a text file with one accession per line (batch mode)",
    )
    parser.add_argument(
        "--max_results",
        type=int,
        default=5,
        help="Max number of papers to fetch per accession (default: 5)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save results as JSON (optional)",
    )
    args = parser.parse_args(argv)

    if args.accession:
        result = run_pipeline(args.accession, max_results=args.max_results)
        if args.output:
            with open(args.output, "w") as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
            print(f"\nResults saved to {args.output}")

    elif args.batch:
        with open(args.batch) as f:
            accessions = [line.strip() for line in f if line.strip()]
        print(f"Batch mode: {len(accessions)} accessions loaded from {args.batch}")
        results = run_batch(accessions, max_results=args.max_results)
        if args.output:
            with open(args.output, "w") as f:
                json.dump(results, f, indent=2, ensure_ascii=False)
            print(f"\nResults saved to {args.output}")
        else:
            print(json.dumps(results, indent=2, ensure_ascii=False))

    else:
        parser.print_help()


if __name__ == "__main__":
    main()