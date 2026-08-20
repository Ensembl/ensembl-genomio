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
"""Vector-search ensemble (BM25 + dense embeddings) over parsed sections, combined with the rule-based result."""

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
from rank_bm25 import BM25Okapi
from sentence_transformers import SentenceTransformer

# Shared ploidy logic — keeps rule-based and vector stages consistent
from .extract import resolve_ploidy_fields, _LEVEL_LABEL

logger = logging.getLogger(__name__)

# ============================================================
# Configuration
# ============================================================

EMBED_MODEL = "pritamdeka/PubMedBERT-mnli-snli-scinli-scitail-mednli-stsb"
INDEX_DIR = "./index_data"
RRF_K = 60

# Reuse section weights from extract.py
SECTION_WEIGHTS = {
    "abstract": 2.0,
    "introduction": 1.5,
    "background": 1.5,
    "results": 1.2,
    "discussion": 1.0,
    "methods": 0.8,
    "materials": 0.8,
    "conclusions": 0.7,
    "supplementary": 0.3,
    "acknowledgements": 0.1,
    "author contributions": 0.0,
    "data availability": 0.0,
}

PRIORITY_SECTIONS = ["abstract", "introduction", "background"]


def get_section_weight(section_name: str) -> float:
    """Return the retrieval weight for a section by matching its name against
    SECTION_WEIGHTS (0.0 = noise section to skip); defaults to 1.0 if unmatched."""
    section_lower = section_name.lower()
    for key, weight in SECTION_WEIGHTS.items():
        if key in section_lower:
            return weight
    return 1.0


# ============================================================
# STEP 1 — Load embedding model
# ============================================================

_embedder: SentenceTransformer | None = None


def get_embedder() -> SentenceTransformer:
    """Return the shared SentenceTransformer, loading it lazily on first use."""
    global _embedder
    if _embedder is None:
        logger.info(f"  [search] Loading embedding model: {EMBED_MODEL}")
        _embedder = SentenceTransformer(EMBED_MODEL)
    return _embedder


# ============================================================
# STEP 2 — Build index from chunks
# ============================================================


def build_index(chunks: list[dict], paper_id: str) -> None:
    """Encode the chunks into dense embeddings and BM25 tokens and persist them
    (embeddings.npy, bm25_tokens.json, metadata.json) under INDEX_DIR."""
    embedder = get_embedder()
    index_dir = Path(INDEX_DIR)
    index_dir.mkdir(parents=True, exist_ok=True)

    # Attach paper_id to each chunk
    for chunk in chunks:
        chunk["paper_id"] = paper_id

    texts = [c["text"] for c in chunks]
    passage_texts = [f"passage: {t}" for t in texts]

    logger.info(f"  [search] Encoding {len(texts)} chunks ...")
    embeddings = embedder.encode(
        passage_texts,
        normalize_embeddings=True,
        show_progress_bar=True,
        batch_size=32,
    )

    # Save dense embeddings
    np.save(index_dir / "embeddings.npy", np.asarray(embeddings).astype("float32"))

    # Save BM25 tokens
    tokenized = [t.lower().split() for t in texts]
    with (index_dir / "bm25_tokens.json").open("w") as handle:
        json.dump(tokenized, handle)

    # Save chunk metadata
    pd.DataFrame(chunks).to_json(
        index_dir / "metadata.json",
        orient="records",
        force_ascii=False,
    )

    logger.info(f"  [search] Index saved to {INDEX_DIR}/ ({len(chunks)} chunks)")


# ============================================================
# STEP 3 — Load index from disk
# ============================================================


def load_index() -> tuple[np.ndarray, BM25Okapi, pd.DataFrame, list[str]]:
    """Load the persisted index from INDEX_DIR and return
    (embeddings, BM25 model, chunk metadata frame, chunk texts)."""
    index_dir = Path(INDEX_DIR)
    embeddings = np.load(index_dir / "embeddings.npy")
    meta = pd.read_json(index_dir / "metadata.json", orient="records")
    with (index_dir / "bm25_tokens.json").open() as handle:
        bm25 = BM25Okapi(json.load(handle))
    texts = meta["text"].tolist()
    logger.info(f"  [search] Index loaded: {len(texts)} chunks")
    return embeddings, bm25, meta, texts


# ============================================================
# STEP 4 — Hybrid search (BM25 + Dense + RRF)
# ============================================================


def rrf_fuse(
    bm25_scores: np.ndarray,
    dense_scores: np.ndarray,
    k: int = 5,
    rrf_k: int = RRF_K,
) -> list[tuple[int, float]]:
    """Fuse BM25 and dense score arrays with Reciprocal Rank Fusion and return the
    top-k (index, fused_score) pairs."""
    bm25_ranks = np.argsort(np.argsort(-bm25_scores))
    dense_ranks = np.argsort(np.argsort(-dense_scores))
    fused = 1.0 / (rrf_k + bm25_ranks) + 1.0 / (rrf_k + dense_ranks)
    top_idx = np.argsort(-fused)[:k]
    return [(int(i), float(fused[i])) for i in top_idx]


def search(
    query: str,
    embeddings: np.ndarray,
    bm25: BM25Okapi,
    meta: pd.DataFrame,
    texts: list[str],
    top_k: int = 5,
) -> list[dict]:
    """Hybrid BM25 + dense search for one query: fuse the two rankings, drop noise
    sections, apply section weights, and return the top_k weighted results."""
    embedder = get_embedder()

    # BM25 scores
    bm25_scores = np.array(bm25.get_scores(query.lower().split()))

    # Dense scores (cosine similarity via dot product — vectors are normalised)
    qvec = embedder.encode([f"query: {query}"], normalize_embeddings=True)[0]
    dense_scores = embeddings @ qvec

    # Hybrid RRF — fetch more candidates before filtering
    hybrid_top = rrf_fuse(bm25_scores, dense_scores, k=top_k * 3)

    results: list = []
    for idx, score in hybrid_top:
        section = meta.iloc[idx]["section"]
        weight = get_section_weight(section)
        if weight == 0.0:
            continue  # skip noise sections entirely
        results.append(
            {
                "rank": len(results) + 1,
                "score": round(score * weight, 4),
                "raw_score": round(score, 4),
                "section": section,
                "text": texts[idx],
                "paper_id": meta.iloc[idx]["paper_id"],
            }
        )

    # Re-rank by weighted score and truncate
    results = sorted(results, key=lambda x: x["score"], reverse=True)[:top_k]
    for i, r in enumerate(results):
        r["rank"] = i + 1
    return results


# ============================================================
# STEP 5 — Infer metadata fields from search results
# ============================================================

SEX_KEYWORDS = ["female", "male", "XY", "ZW", "pistillate", "staminate", "dioecious"]

# Evidence weighting — must match extract.py so both stages agree
_COMPOUND_BOOST = 6.0
_FORMULA_BOOST = 10.0


def infer_ploidy(results: list) -> dict:
    """Ploidy (INSDC schema tuple) from retrieved chunks, using the SAME
    resolver as the rule-based stage (extract.resolve_ploidy_fields) so both
    stages apply identical ancestor-aware, compound-aware logic.
    """
    segments = []
    for r in results:
        w = get_section_weight(r["section"])
        if w == 0.0:
            continue
        is_priority = any(p in r["section"].lower() for p in PRIORITY_SECTIONS)
        segments.append((r["text"], w, is_priority))
    return resolve_ploidy_fields(segments)


def infer_sex(results: list[dict]) -> dict:
    """Infer the sequenced individual's sex from retrieved chunks via SEX_KEYWORDS,
    returning the match with context or {"value": "unknown"}."""
    for r in results:
        if get_section_weight(r["section"]) == 0.0:
            continue
        text_lower = r["text"].lower()
        for kw in SEX_KEYWORDS:
            if re.search(r"\b" + re.escape(kw.lower()) + r"\b", text_lower):
                idx = text_lower.find(kw.lower())
                start = max(0, idx - 60)
                end = min(len(r["text"]), idx + 60)
                return {
                    "value": kw,
                    "section": r["section"],
                    "context": r["text"][start:end],
                }
    return {"value": "unknown", "evidence": None}


def infer_species(results: list[dict]) -> dict:
    """Return the best species-bearing chunk, preferring a priority section
    (abstract/introduction/background) then the top-scored result."""
    for r in results:
        if any(p in r["section"].lower() for p in PRIORITY_SECTIONS):
            return {
                "text": r["text"],
                "section": r["section"],
                "score": r["score"],
            }
    if results:
        return {
            "text": results[0]["text"],
            "section": results[0]["section"],
            "score": results[0]["score"],
        }
    return {"text": None, "section": None}


def infer_strain(results: list[dict]) -> dict:
    """Return the top retrieved chunk likely to name a strain/cultivar, skipping
    non-informative sections (author contributions, conclusions, acknowledgements)."""
    EXCLUDE = ["author contribution", "conclusions", "acknowledgement"]
    for r in results:
        if any(e in r["section"].lower() for e in EXCLUDE):
            continue
        return {
            "text": r["text"],
            "section": r["section"],
            "score": r["score"],
        }
    if results:
        return {"text": results[0]["text"], "section": results[0]["section"]}
    return {"text": None, "section": None}


# ============================================================
# STEP 6 — Ensemble: combine rule-based + vector results
# ============================================================

_MECH_SPEC = {"amphidiploid": 3, "segmental_allopolyploid": 3, "allopolyploid": 1, "autopolyploid": 1}


def _merge_mechanism(a: str | None, b: str | None) -> str | None:
    cands = [m for m in (a, b) if m]
    if not cands:
        return None
    return max(cands, key=lambda m: _MECH_SPEC.get(m, 0))


def ensemble_ploidy(rule_result: dict, vector_result: dict) -> dict:
    """Combine rule-based and vector ploidy on the `level` field, merging the
    other schema elements (mechanism / irregular / derivation)."""
    r_lvl = rule_result.get("level")
    v_lvl = vector_result.get("level")
    r_conf = rule_result.get("confidence", 0.0)
    v_conf = vector_result.get("confidence", 0.0)

    mechanism = _merge_mechanism(rule_result.get("mechanism"), vector_result.get("mechanism"))
    irregular = rule_result.get("irregular") or vector_result.get("irregular") or False
    derivation = rule_result.get("derivation") or vector_result.get("derivation")

    def out(
        level: int | None, conf: float, method: str, evidence: object, origin: dict | None = None
    ) -> dict:
        mech = mechanism
        if level == 1:
            mech = None
        elif level is not None and level < 3:
            mech = None
        # carry the status/source from whichever result supplied the level
        if origin is not None and (origin.get("status") or origin.get("source")):
            status = origin.get("status")
            source = origin.get("source")
        elif level is not None or mech is not None:
            status, source = "found_in_paper", "paper"
        else:
            status, source = "not_stated_in_paper", None
        return {
            "level": level,
            "mechanism": mech,
            "irregular": irregular,
            "derivation": derivation,
            "tuple": [level, mech, irregular, derivation],
            "level_label": _LEVEL_LABEL.get(level) if level is not None else None,
            "status": status,
            "source": source,
            "confidence": conf,
            "method": method,
            "evidence": evidence,
        }

    if r_lvl is None and v_lvl is None:
        if mechanism:
            return out(
                None,
                max(r_conf, v_conf, 0.4),
                "mechanism_only",
                rule_result.get("evidence") or vector_result.get("evidence"),
                rule_result,
            )
        return out(None, 0.0, "both_unknown", None, rule_result)
    elif r_lvl == v_lvl:
        return out(
            r_lvl,
            min(round((r_conf + v_conf) / 2 + 0.1, 2), 0.95),
            "consensus",
            rule_result.get("evidence"),
            rule_result,
        )
    elif r_lvl is None:
        return out(v_lvl, round(v_conf * 0.9, 2), "vector_only", vector_result.get("evidence"), vector_result)
    elif v_lvl is None:
        return out(r_lvl, round(r_conf * 0.9, 2), "rule_only", rule_result.get("evidence"), rule_result)
    else:
        # Conflict — penalise confidence, pick higher-confidence level
        if r_conf >= v_conf:
            return out(
                r_lvl,
                round(r_conf * 0.8, 2),
                "rule_wins (conflict)",
                rule_result.get("evidence"),
                rule_result,
            )
        else:
            return out(
                v_lvl,
                round(v_conf * 0.8, 2),
                "vector_wins (conflict)",
                vector_result.get("evidence"),
                vector_result,
            )


def ensemble_cultivars(
    rule_cultivars: list,
    vector_strain: dict,
    sections: dict,
) -> dict:
    """Split rule-based cultivar candidates into confirmed vs unconfirmed by
    checking each against the full section text; returns both lists plus section."""
    full_text = " ".join(text for sec, text in sections.items() if get_section_weight(sec) > 0)
    confirmed, unconfirmed = [], []
    for c in rule_cultivars:
        if len(c) > 30:
            continue
        if re.search(re.escape(c.lower()), full_text.lower()):
            confirmed.append(c)
        else:
            unconfirmed.append(c)
    return {
        "confirmed": confirmed,
        "unconfirmed": unconfirmed,
        "section": vector_strain.get("section"),
    }


def ensemble_sex(rule_sex: dict, vector_sex: dict) -> dict:
    """Combine rule-based and vector sex calls: consensus when they agree, the
    non-unknown one when only one has a value, else 'unknown' (conflict)."""
    r_val = rule_sex.get("value", "unknown")
    v_val = vector_sex.get("value", "unknown")
    if r_val == v_val:
        return {"value": r_val, "method": "consensus"}
    elif r_val == "unknown":
        return {"value": v_val, "method": "vector_only"}
    elif v_val == "unknown":
        return {"value": r_val, "method": "rule_only"}
    else:
        return {"value": "unknown", "method": "conflict"}


# ============================================================
# MAIN FUNCTION — run all steps above in order
#           input:  parsed dict from parse.py
#                   rule_result from extract.py
#           output: final ensemble metadata
# ============================================================

QUERIES = {
    "ploidy": "ploidy level diploid triploid tetraploid hexaploid octoploid "
    "polyploid allopolyploid autopolyploid amphidiploid "
    "whole genome duplication WGD karyotype "
    "chromosome number chromosome count 2n 4x 6x",
    "strain": "cultivar strain ecotype variety accession",
    "sex": "female male sex chromosome XY ZW dioecious",
    "species": "scientific name genus species organism",
}


def run_vector_search(parsed: dict, rule_result: dict) -> dict:
    """Build the per-paper index, run hybrid search for each metadata field, infer
    ploidy/sex/species/strain from the results, and ensemble them with the
    rule-based result into the final metadata dict."""
    sections = parsed.get("sections", {})
    chunks = parsed.get("chunks", [])
    paper_id = parsed.get("source", "unknown")

    # Build index for this paper
    logger.info("\n  [search] Building index ...")
    build_index(chunks, paper_id=paper_id)

    # Load index
    embeddings, bm25, meta, texts = load_index()

    # Search for each metadata field
    logger.info("\n  [search] Running hybrid search ...")
    vector_results = {}
    for field, query in QUERIES.items():
        results = search(query, embeddings, bm25, meta, texts, top_k=5)
        vector_results[field] = results

    # Infer from search results
    v_ploidy = infer_ploidy(vector_results["ploidy"])
    v_sex = infer_sex(vector_results["sex"])
    v_species = infer_species(vector_results["species"])
    v_strain = infer_strain(vector_results["strain"])

    # Ensemble with rule-based results
    logger.info("\n  [search] Running ensemble ...")
    final: dict = {
        # NCBI-sourced fields (pass through from extract.py)
        "assembly_accession": rule_result.get("assembly_accession"),
        "assembly_name": rule_result.get("assembly_name"),
        "taxon_id": rule_result.get("taxon_id"),
        "chromosome_number": rule_result.get("chromosome_number"),
        # Ensemble fields
        "ploidy": ensemble_ploidy(rule_result["ploidy"], v_ploidy),
        "cultivars": ensemble_cultivars(rule_result["cultivars"], v_strain, sections),
        "sex": ensemble_sex(rule_result["sex"], v_sex),
        # Species: NCBI > rule-based > vector
        "species": rule_result.get("species")
        or {
            "value": v_species.get("text", "unknown"),
            "source": "vector_search",
        },
    }

    # Retain the ploidy evidence passages (semantic-retrieval provenance):
    # the top ploidy-query chunks, each with its section and fused score, so a
    # curator can see exactly which passages support the ploidy call.
    final["ploidy"]["evidence_passages"] = [
        {
            "text": r["text"],
            "section": r["section"],
            "score": r["score"],
        }
        for r in vector_results["ploidy"][:5]
    ]

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("Final Ensemble Results")
    logger.info("=" * 60)
    logger.info(f"  assembly_accession = {final['assembly_accession']}")
    logger.info(f"  assembly_name      = {final['assembly_name']}")
    logger.info(f"  taxon_id           = {final['taxon_id']}")
    logger.info(
        f"  species            = {final['species']['value']} " f"(source: {final['species']['source']})"
    )
    p = final["ploidy"]
    logger.info(
        f"  ploidy             = {p['tuple']} "
        f"(level={p['level']}, confidence: {p['confidence']}, method: {p['method']})"
    )
    c = final["chromosome_number"]
    logger.info(f"  chromosome_number  = {c['value']} (source: {c['source']})")
    logger.info(f"  cultivars          = {final['cultivars']['confirmed']}")
    logger.info(f"  sex                = {final['sex']['value']} " f"(method: {final['sex']['method']})")
    logger.info("=" * 60)
    logger.info("\n  Confidence guide:")
    logger.info("    consensus         both methods agree     -> highest reliability")
    logger.info("    rule_only         rule-based only        -> needs review")
    logger.info("    vector_only       vector search only     -> needs review")
    logger.info("    conflict          methods disagree       -> manual check required")

    return final
