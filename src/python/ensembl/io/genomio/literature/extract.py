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
"""Rule-based extraction of ploidy, chromosome number, cultivar/strain and sex from parsed text."""

import logging
import re
from typing import Optional

logger = logging.getLogger(__name__)

# ============================================================
# Section weights
# ============================================================

SECTION_WEIGHTS = {
    "title": 2.5,  # paper title — most authoritative for ploidy
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
    "author contributions": 0.0,  # noise — skip entirely
    "data availability": 0.0,  # noise — skip entirely
}


def get_section_weight(section_name: str) -> float:
    """Return the ploidy-evidence weight for a section by matching its name against
    SECTION_WEIGHTS (0.0 = noise section to skip); defaults to 1.0 if unmatched."""
    section_lower = section_name.lower()
    for key, weight in SECTION_WEIGHTS.items():
        if key in section_lower:
            return weight
    return 1.0  # default


PRIORITY_SECTIONS = ["abstract", "introduction", "background"]


# ============================================================
# STEP 1 — Extract ploidy  (INSDC Ploidy Schema v1.1)
# ============================================================
#
# ploidy = (level, mechanism, irregular, derivation)
#
#   level       int | None     number of complete haploid genome sets (x)
#   mechanism   str | None      origin of polyploidy (controlled vocab)
#   irregular   False | str     deviation from a clean ploidy state
#   derivation  str | None      biotechnological derivation of the line
#
# Design notes:
#   - `level` is the x-based integer (diploid = 2), NOT the gametic number n.
#   - The paper title is the most authoritative source and is injected as a
#     high-weight section by extract_metadata().
#   - Polyploid papers mention their DIPLOID ANCESTORS more often than the
#     organism's own ploidy, so mentions in an ancestor/progenitor/donor/
#     "other diploid" context are excluded from the vote.
#   - "doubled haploid" etc. are NOT discarded: they populate `derivation`.
#   - "amphidiploid"/"dihaploid" are excluded from the level word match
#     (their embedded "diploid"/"haploid" would give a wrong level) and are
#     handled by the mechanism / derivation scanners instead.
# ============================================================

# ---- level: ploidy word -> integer ----
_WORD_LEVEL = {
    "monoploid": 1,
    "diploid": 2,
    "triploid": 3,
    "tetraploid": 4,
    "pentaploid": 5,
    "hexaploid": 6,
    "heptaploid": 7,
    "octoploid": 8,
    "decaploid": 10,
    "dodecaploid": 12,
}
_LEVEL_LABEL = {v: k for k, v in _WORD_LEVEL.items()}

# Level words with optional allo/auto prefix. Lookbehind/lookahead exclude
# "amphidiploid" (preceded by "amphi"), "diploidization", and require a real
# word boundary; trailing y/s allow "diploids" / "diploidy". "haploid" is
# deliberately NOT a level word (ambiguous: gametic phase / assembly type).
_LEVEL_RE = re.compile(
    r"(?<![A-Za-z])(allo|auto)?[\-\s]?"
    r"(monoploid|diploid|triploid|tetraploid|pentaploid|hexaploid|heptaploid|octoploid|decaploid|dodecaploid)"
    r"(?:y|s)?(?![A-Za-z])",
    re.IGNORECASE,
)

# 2n = Nx = M  ->  N is the ploidy level (most reliable)
_CHROM_FORMULA_RE = re.compile(r"2n\s*=\s*(\d+)\s*x\s*=\s*\d+", re.IGNORECASE)

# Contexts where a diploid/polyploid term describes an ancestor/relative,
# not the assembled organism -> excluded from the organism-level vote.
_ANCESTOR_CONTEXT_RE = re.compile(
    r"ancestor|ancestral|progenitor|parental|\brelativ|\bdonor|"
    r"sub-?genome|derived\s+from|descend|two\s+diploid|both\s+diploid|"
    r"other\s+diploid|related\s+diploid|diploid\s+relativ|"
    r"diploid\s+ancestor|diploid\s+progenitor|diploid\s+wild|wild\s+diploid|"
    r"each\s+diploid",
    re.IGNORECASE,
)

# ---- mechanism (Element 2): (name, pattern, specificity_rank) ----
_MECH_PATTERNS = [
    ("amphidiploid", re.compile(r"amphidiploid", re.IGNORECASE), 3),
    ("segmental_allopolyploid", re.compile(r"segmental[\s_\-]+allo\w*ploid", re.IGNORECASE), 3),
    ("allopolyploid", re.compile(r"\ballo[\-\s]?\w*ploid", re.IGNORECASE), 1),
    ("autopolyploid", re.compile(r"\bauto[\-\s]?\w*ploid", re.IGNORECASE), 1),
]

# ---- irregular (Element 3) ----
_IRREGULAR_PATTERNS = [
    ("aneuploid", re.compile(r"aneuploid", re.IGNORECASE)),
    ("mixoploid", re.compile(r"mixoploid", re.IGNORECASE)),
    ("endopolyploid", re.compile(r"endopolyploid|endoreduplicat", re.IGNORECASE)),
]

# ---- derivation (Element 4) ----
_DERIVATION_PATTERNS = [
    ("doubled_haploid", re.compile(r"double[d]?[\s\-]+haploid", re.IGNORECASE)),
    ("dihaploid", re.compile(r"dihaploid", re.IGNORECASE)),
    (
        "induced_polyploid",
        re.compile(
            r"induced[\s\-]+polyploid|colchicine[\s\-]+(?:doubl|induc|treat)|"
            r"artificial(?:ly)?[\s\-]+induced[\s\-]+polyploid",
            re.IGNORECASE,
        ),
    ),
]

# Evidence weighting
_COMPOUND_BOOST = 6.0
_FORMULA_BOOST = 10.0


def find_level_mentions(text: str, window: int = 50) -> list:
    """Find every ploidy-LEVEL mention and classify it.

    Shared by extract.py and search.py. Each item:
        {level, prefix, is_compound, is_ancestor, is_formula, context}
    """
    mentions = []

    # 2n = Nx = M  (definitive)
    for m in _CHROM_FORMULA_RE.finditer(text):
        x = int(m.group(1))
        if not (1 <= x <= 20):
            continue
        ctx = text[max(0, m.start() - window) : m.end() + window]
        mentions.append(
            {
                "level": x,
                "prefix": None,
                "is_compound": True,
                "is_ancestor": False,
                "is_formula": True,
                "context": ctx,
            }
        )

    # ploidy words
    for m in _LEVEL_RE.finditer(text):
        prefix = (m.group(1) or "").lower()
        base = m.group(2).lower()
        level = _WORD_LEVEL[base]
        ctx = text[max(0, m.start() - window) : m.end() + window]
        is_compound = prefix in ("allo", "auto")
        is_ancestor = (not is_compound) and bool(_ANCESTOR_CONTEXT_RE.search(ctx))
        mentions.append(
            {
                "level": level,
                "prefix": prefix,
                "is_compound": is_compound,
                "is_ancestor": is_ancestor,
                "is_formula": False,
                "context": ctx,
            }
        )

    return mentions


# Mechanism words almost always describe the focal organism. Suppress only
# when an ancestor noun sits IMMEDIATELY next to the term (e.g. "allopolyploid
# ancestor"), not merely somewhere in the sentence ("allohexaploid arising
# from a tetraploid progenitor" still describes the organism).
_MECH_ANCESTOR_RE = re.compile(r"ancestor|ancestral|progenitor|\brelativ|\bdonor", re.IGNORECASE)


def detect_mechanism(text: str, window: int = 22) -> Optional[str]:
    """Return the most specific organism-level polyploidy mechanism, or None."""
    best, best_rank = None, -1
    for name, rx, rank in _MECH_PATTERNS:
        m = rx.search(text)
        if not m:
            continue
        ctx = text[max(0, m.start() - window) : m.end() + window]
        if _MECH_ANCESTOR_RE.search(ctx):
            continue  # the term itself qualifies an ancestor, not the organism
        if rank > best_rank:
            best, best_rank = name, rank
    return best


def detect_irregular(text: str) -> str | bool:
    """Return the ploidy-irregularity term (aneuploid / mixoploid / endopolyploid)
    found in the text, or False if none is present."""
    for name, rx in _IRREGULAR_PATTERNS:
        if rx.search(text):
            return name
    return False


def detect_derivation(text: str) -> Optional[str]:
    """Return the biotechnological derivation of the line (doubled_haploid /
    dihaploid / induced_polyploid) found in the text, or None."""
    for name, rx in _DERIVATION_PATTERNS:
        if rx.search(text):
            return name
    return None


def resolve_ploidy_fields(segments: list) -> dict:
    """Resolve the 4-element ploidy tuple from weighted text segments.

    segments: iterable of (text, weight, is_priority).
    Returns the schema dict (level / mechanism / irregular / derivation + meta).
    Used by both extract.py (sections) and search.py (retrieved chunks).
    """
    level_counts: dict = {}
    best_ev: dict = {}
    mech_votes: dict = {}
    irregular: str | bool = False
    derivation = None
    n_org = 0

    for text, weight, is_priority in segments:
        if not text or weight == 0.0:
            continue

        for men in find_level_mentions(text):
            if men["is_ancestor"]:
                continue
            n_org += 1
            if men["is_formula"]:
                w = weight * _FORMULA_BOOST
            elif men["is_compound"]:
                w = weight * _COMPOUND_BOOST
            else:
                w = weight
            lvl = men["level"]
            level_counts[lvl] = level_counts.get(lvl, 0.0) + w

            rank_new = (men["is_formula"], men["is_compound"], is_priority)
            cur = best_ev.get(lvl)
            if cur is None or rank_new > cur["_rank"]:
                best_ev[lvl] = {
                    "context": men["context"],
                    "is_priority": is_priority,
                    "is_formula": men["is_formula"],
                    "_rank": rank_new,
                }

        m = detect_mechanism(text)
        if m:
            mech_votes[m] = mech_votes.get(m, 0.0) + weight
        if not irregular:
            irregular = detect_irregular(text) or False
        if derivation is None:
            derivation = detect_derivation(text)

    # ----- choose level -----
    if level_counts:
        level = max(level_counts, key=lambda k: level_counts[k])
        total = sum(level_counts.values())
        base_conf = level_counts[level] / total
        ev = best_ev.get(level, {})
        priority_bonus = 0.1 if ev.get("is_priority") else 0.0
        confidence = min(round(base_conf + priority_bonus, 2), 0.95)
        method = "chromosome_formula" if ev.get("is_formula") else "text_weighted_vote"
        evidence = ev.get("context", "")
    else:
        level = None
        confidence = 0.0
        method = "no_text_evidence"
        evidence = ""

    # ----- choose mechanism -----
    mechanism = None
    if mech_votes:
        spec = {"amphidiploid": 3, "segmental_allopolyploid": 3, "allopolyploid": 1, "autopolyploid": 1}
        mechanism = max(mech_votes, key=lambda k: (spec.get(k, 0), mech_votes[k]))
        if confidence == 0.0:
            confidence = 0.5  # polyploid confirmed but level not stated

    # ----- schema validation / normalisation -----
    if level == 1:
        mechanism = None  # monoploid MUST have null mechanism
    elif level is not None and level < 3:
        mechanism = None  # diploid SHOULD have null mechanism

    if level is not None or mechanism is not None:
        status, source = "found_in_paper", "paper"
    else:
        status, source = "not_stated_in_paper", None

    return {
        "level": level,
        "mechanism": mechanism,
        "irregular": irregular,
        "derivation": derivation,
        "tuple": [level, mechanism, irregular, derivation],
        "level_label": _LEVEL_LABEL.get(level) if level is not None else None,
        "status": status,
        "source": source,
        "confidence": confidence,
        "method": method,
        "evidence": evidence,
        "evidence_count": n_org,
        "weighted_counts": {k: round(v, 2) for k, v in level_counts.items()},
    }


def extract_ploidy(weighted_sections: dict) -> dict:
    """Extract the INSDC ploidy tuple from weighted sections by delegating to
    resolve_ploidy_fields (marking abstract/introduction/background as priority)."""
    segments = [
        (text, weight, any(p in sec.lower() for p in PRIORITY_SECTIONS))
        for sec, (text, weight) in weighted_sections.items()
    ]
    return resolve_ploidy_fields(segments)


# ============================================================
# Reference ploidy fallback (OPTIONAL, transparently sourced)
# ============================================================
# Used ONLY when neither the paper text nor NCBI yields a level.
# Returns the schema tuple with method="reference_taxon" and modest
# confidence so it is visibly NOT derived from the paper.
#
# OPEN QUESTION (pending mentor): the proposal's anti-hallucination rule
# argues for level=None when the text is silent. To switch to that behaviour,
# delete KNOWN_PLOIDY and the fallback block in extract_metadata().
# ============================================================

KNOWN_PLOIDY = {
    "3702": 2,  # Arabidopsis thaliana
    "4081": 2,  # Solanum lycopersicum (tomato)
    "4577": 2,  # Zea mays
    "39947": 2,  # Oryza sativa Japonica
    "3818": 4,  # Arachis hypogaea (peanut)
    "4565": 6,  # Triticum aestivum (bread wheat)
}


# ============================================================
# STEP 2 — Extract chromosome number from text
# ============================================================

# Order matters: the compound "2n = Nx = M" form must be tried BEFORE the bare
# "2n = N" form, otherwise "2n = 4x = 52" matches "2n = 4" and captures the ploidy
# coefficient (4) instead of the chromosome count (52).
CHROMOSOME_PATTERNS = [
    r"2n\s*=\s*\d+x\s*=\s*(\d+)",  # 2n = 4x = 52  -> 52 (total count)
    r"2n\s*=\s*(\d+)",  # 2n = 48
    r"n\s*=\s*(\d+)",  # n = 24 (haploid)
    r"(\d+)\s+chromosome",  # 48 chromosomes
    r"chromosome\s+number\s+(?:of\s+)?(\d+)",
]


def extract_chromosome_number(weighted_sections: dict) -> dict:
    """Extract a chromosome number from text, scanning sections highest-weight
    first and returning the first CHROMOSOME_PATTERNS match with its context."""
    priority_order = sorted(
        weighted_sections.items(), key=lambda x: x[1][1], reverse=True  # sort by weight descending
    )

    for sec_name, (text, weight) in priority_order:
        if weight == 0.0:
            continue
        for pattern in CHROMOSOME_PATTERNS:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                number = match.group(1)
                ctx_start = max(0, match.start() - 50)
                ctx_end = min(len(text), match.end() + 50)
                return {
                    "value": number,
                    "evidence": text[ctx_start:ctx_end],
                    "section": sec_name,
                    "source": "text_extraction",
                }

    return {"value": None, "evidence": "", "source": "text_extraction"}


# ============================================================
# STEP 3 — Extract cultivar / strain names from text
# ============================================================

CULTIVAR_PATTERNS = [
    r"(?i:cultivar)\s+['\u2018\u2019\"]([^'\u2018\u2019\"]+)['\u2018\u2019\"]",
    r"(?i:cv)\.\s*['\u2018\u2019\"]?([A-Z][a-zA-Z\s]+)",
    r"['\u2018\u2019]([A-Z][a-zA-Z\s]{2,25})['\u2018\u2019]",
    r"(?i:cultivar)[:\s]+([A-Z][\w\-]*(?: [A-Z0-9][\w\-]*){0,2})",
    r"(?i:strain)[:\s]+([A-Z][\w\-]*(?: [A-Z0-9][\w\-]*){0,2})",
    r"(?i:ecotype)[:\s]+([A-Z][\w\-]*(?: [A-Z0-9][\w\-]*){0,2})",
]

CULTIVAR_STOPWORDS = {
    "the",
    "and",
    "for",
    "this",
    "that",
    "with",
    "from",
    "table",
    "figure",
    "supplementary",
    "data",
    "note",
    "results",
    "methods",
    "background",
    "abstract",
    "asm",
    "bam",
    "fasta",
    "fastq",
    "csv",
    "tsv",
    "json",
    "xml",
    "hifiasm",
    "busco",
    "samtools",
    "orthofinder",
    "blast",
    "minimap",
    "viridiplantae",
    "embryophyta",
    "rosaceae",
    "eudicots",
    "haplome",
    "haplotype",
    "genbank",
    "refseq",
    "bioproject",
}


def _looks_like_name(name: str) -> bool:
    if not (2 < len(name) <= 30):
        return False
    tokens = name.split()
    if not (1 <= len(tokens) <= 3):
        return False
    return all(token[0].isupper() or token[0].isdigit() for token in tokens)


def extract_cultivars(filtered_text: str) -> list:
    """Extract candidate cultivar/strain/ecotype names via CULTIVAR_PATTERNS,
    dropping stopwords and non-name-like tokens; returns a sorted unique list."""
    cultivars = set()
    for pattern in CULTIVAR_PATTERNS:
        for match in re.findall(pattern, filtered_text):
            name = match.strip().rstrip(".").strip()
            if name.lower() in CULTIVAR_STOPWORDS:
                continue
            if _looks_like_name(name):
                cultivars.add(name)
    return sorted(cultivars)


# ============================================================
# STEP 4 — Extract sex from text
# ============================================================

SEX_PATTERNS = [
    r"(female|male)\s+(individual|specimen|plant|animal|organism)",
    r"sex[:\s]+(male|female|hermaphrodite)",
    r"(maternal|paternal)\s+",
]


def extract_sex(weighted_sections: dict) -> dict:
    """Extract the sequenced individual's sex via SEX_PATTERNS, scanning non-noise
    sections; returns {"value": "unknown"} if not stated."""
    for sec_name, (text, weight) in weighted_sections.items():
        if weight == 0.0:
            continue
        for pattern in SEX_PATTERNS:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                return {
                    "value": match.group(1).lower(),
                    "section": sec_name,
                    "evidence": text[max(0, match.start() - 30) : match.end() + 30],
                }
    return {"value": "unknown", "evidence": ""}


# ============================================================
# STEP 5 — Extract species name from text
# ============================================================


def extract_species(sections: dict, weighted_sections: dict) -> Optional[str]:
    """Extract a binomial species name from text, preferring the abstract then
    high-weight sections; returns "unknown" if no binomial is found. Used only as
    a fallback when NCBI provides no scientific_name."""
    pattern1 = r"([A-Z][a-z]+\s+[a-z]+)(?:\s+(?:Borkh\.|L\.|var\.))"
    pattern2 = r"([A-Z][a-z]{2,15}\s+[a-z]{2,15})"

    abstract = sections.get("Abstract", "")
    for pat in [pattern1, pattern2]:
        m = re.search(pat, abstract)
        if m:
            return m.group(1)

    for sec_name, (text, weight) in weighted_sections.items():
        if weight < 1.5:
            continue
        for pat in [pattern1, pattern2]:
            m = re.search(pat, text)
            if m:
                return m.group(1)

    filtered = " ".join(t for _, (t, w) in weighted_sections.items() if w > 0)
    m = re.search(pattern2, filtered)
    if m:
        return m.group(1)

    return "unknown"


# ============================================================
# MAIN FUNCTION — run all steps above in order
#           input:  parsed dict from parse.py
#                   assembly dict from fetch.py
#                   paper_title (optional, high-signal for ploidy)
#           output: all extracted metadata fields
# ============================================================


def extract_metadata(parsed: dict, assembly: dict, paper_title: Optional[str] = None) -> dict:
    """Run the full rule-based extraction for one paper: ploidy, chromosome number,
    cultivars, sex and species. NCBI assembly fields take precedence; text is used
    to fill gaps, and a clearly-sourced reference ploidy is a last resort when the
    paper is silent. Returns the combined metadata dict."""
    sections = dict(parsed.get("sections", {}))

    # Paper title is the most authoritative source for ploidy and is often
    # NOT part of the fetched abstract/fulltext body, so inject it.
    if paper_title:
        sections["Title"] = paper_title

    weighted_sections = {sec: (text, get_section_weight(sec)) for sec, text in sections.items()}

    filtered_text = " ".join(text for sec, text in sections.items() if get_section_weight(sec) > 0)

    logger.info("  [extract] Running rule-based extraction ...")

    # ── Fields from NCBI ──────────────────────────────────────
    ncbi_accession = assembly.get("assembly_accession")
    ncbi_name = assembly.get("assembly_name")
    ncbi_taxon = assembly.get("taxon_id")
    ncbi_species = assembly.get("scientific_name")
    ncbi_chrom = assembly.get("chromosome_number")
    ncbi_chrom_source = assembly.get("chromosome_source")

    if ncbi_species:
        species = {"value": ncbi_species, "source": "ncbi"}
    else:
        species = {
            "value": extract_species(sections, weighted_sections),
            "source": "text_extraction",
        }

    if ncbi_chrom:
        chromosome_number = {"value": str(ncbi_chrom), "source": ncbi_chrom_source}
    else:
        logger.info("  [extract] chromosome_number not in NCBI → extracting from text ...")
        chromosome_number = extract_chromosome_number(weighted_sections)

    # ── Text-only fields ──────────────────────────────────────
    ploidy = extract_ploidy(weighted_sections)
    cultivars = extract_cultivars(filtered_text)
    sex = extract_sex(weighted_sections)

    # ── Ploidy reference fallback (paper is silent on ploidy) ────────────
    #     Mentor: a reference value is acceptable IF it is clearly sourced and
    #     kept separate from paper-derived values. Prefer GoaT (taxon-based,
    #     attached by fetch.enrich_*), fall back to the small built-in table.
    if ploidy["level"] is None:
        ref_level = assembly.get("reference_ploidy")
        ref_source = "GoaT" if ref_level is not None else None
        if ref_level is None and str(ncbi_taxon) in KNOWN_PLOIDY:
            ref_level, ref_source = KNOWN_PLOIDY[str(ncbi_taxon)], "reference_db"
        if ref_level is not None:
            logger.info(
                f"  [extract] ploidy not stated in paper → reference: "
                f"{ref_level} (taxon {ncbi_taxon}, source {ref_source})"
            )
            ploidy = {
                "level": ref_level,
                "mechanism": None,
                "irregular": False,
                "derivation": None,
                "tuple": [ref_level, None, False, None],
                "level_label": _LEVEL_LABEL.get(ref_level),
                "status": "from_reference",
                "source": ref_source,
                "confidence": 0.6,
                "method": "reference_taxon",
                "evidence": f"(reference {ref_source}: taxon {ncbi_taxon} documented as level {ref_level})",
                "evidence_count": 0,
                "weighted_counts": {},
            }
        else:
            logger.info(
                f"  [extract] ploidy not stated in paper and no reference available "
                f"→ status: not_stated_in_paper"
            )

    logger.info(f"  [extract] Done.")
    logger.info(f"            species           = {species['value']} (source: {species['source']})")
    logger.info(
        f"            ploidy            = {ploidy['tuple']} "
        f"(level={ploidy['level']}, confidence: {ploidy['confidence']}, method: {ploidy['method']})"
    )
    logger.info(
        f"            chromosome_number = {chromosome_number['value']} (source: {chromosome_number['source']})"
    )
    logger.info(f"            cultivars         = {cultivars}")
    logger.info(f"            sex               = {sex['value']}")

    return {
        "assembly_accession": ncbi_accession,
        "assembly_name": ncbi_name,
        "taxon_id": ncbi_taxon,
        "species": species,
        "ploidy": ploidy,
        "chromosome_number": chromosome_number,
        "cultivars": cultivars,
        "sex": sex,
    }
