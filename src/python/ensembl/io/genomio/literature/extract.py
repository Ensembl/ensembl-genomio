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
_CHROM_FORMULA_RE = re.compile(r"2n\s*=\s*(\d+)\s*x\s*=\s*(\d+)", re.IGNORECASE)

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

# Technical uses of a level word that say nothing about the organism's ploidy:
# "diploid assembly" / "diploid-aware assembler" (assembly method) and "diploid
# chromosome number" (= 2n, used for polyploids too). Matched right after the term.
_TECHNICAL_CONTEXT_RE = re.compile(
    r"^[\s\-]*(?:genome\s+)?(?:assembl|chromosome\s+number|number\b|mode\b|aware|phas)",
    re.IGNORECASE,
)

# Species names near a level word. Abbreviated binomials ("D. nipponica") are
# unambiguous; full binomials need a Latin-looking epithet and a genus that is not
# a common capitalised word, to avoid matching prose such as "Our data".
_ABBREV_SPECIES_RE = re.compile(r"\b([A-Z])\.\s?([a-z]{3,})\b")
_FULL_SPECIES_RE = re.compile(
    r"\b([A-Z][a-z]{2,})\s([a-z]{3,}(?:a|um|us|is|ii|ae|i|ensis|oides|ens|ana|ica|ata))\b"
)
_NOT_GENUS = {
    "the",
    "this",
    "these",
    "that",
    "those",
    "our",
    "here",
    "for",
    "from",
    "with",
    "each",
    "both",
    "all",
    "some",
    "such",
    "most",
    "many",
    "table",
    "figure",
    "fig",
    "data",
    "note",
    "section",
    "supplementary",
    "they",
    "its",
    "their",
    "when",
    "while",
    "after",
    "before",
}
_NOT_EPITHET = {"data", "area", "formula", "media", "criteria", "phenomena", "schema", "via", "idea", "extra"}


def target_names(scientific_name: str | None, common_name: str | None = None) -> list[str]:
    """Lower-cased names that identify the assembled organism in text: genus,
    specific epithet and the head noun of the common name (e.g. ["triticum",
    "aestivum", "wheat"]). Only the head noun is used: modifiers such as "common"
    or "domestic" would match unrelated text."""
    words = (scientific_name or "").split()[:2]
    head = re.split(r"[\s\-]+", (common_name or "").strip())[-1]
    if len(head) > 2:
        words.append(head)
    return [w.lower() for w in words if w]


def _names_other_species(context: str, names: list[str]) -> bool:
    """True if `context` names a species but none of the target organism's names,
    i.e. a nearby ploidy statement most likely describes another organism."""
    if not names:
        return False
    lowered = context.lower()
    if any(re.search(rf"\b{re.escape(name)}", lowered) for name in names):
        return False
    genus = names[0]
    for initial, epithet in _ABBREV_SPECIES_RE.findall(context):
        if initial.lower() != genus[:1] or epithet.lower() not in names:
            return True
    for genus_word, epithet in _FULL_SPECIES_RE.findall(context):
        if genus_word.lower() not in _NOT_GENUS and epithet.lower() not in _NOT_EPITHET:
            return True
    return False


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
    ("doubled_monoploid", re.compile(r"double[d]?[\s\-]+monoploid", re.IGNORECASE)),
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

# Organism-level "haploid" ("a haploid yeast", "haploid and diploid generations").
# "haploid" is not a reportable level (it is ambiguous), but an organism-level
# mention is counter-evidence against reporting diploid. Technical uses ("haploid
# genome size", "haploid chromosome number", "haploid assembly") and doubled /
# pseudo haploids are ignored.
_HAPLOID_RE = re.compile(r"(?<![A-Za-z\-])haploid(?:s|y)?(?![A-Za-z\-])", re.IGNORECASE)
_HAPLOID_SKIP_BEFORE_RE = re.compile(r"(?:double[d]?|pseudo)[\s\-]*$", re.IGNORECASE)
_HAPLOID_TECHNICAL_RE = re.compile(
    r"^[\s\-]*(?:genome|chromosome|number|assembl|set|size|content|dna|c-value|sequence|version|"
    r"consensus|represent|resolved|phas|mode|reference|length|coverage)",
    re.IGNORECASE,
)

# A supporting mention is "anchored" when a target name occurs this close to it,
# or when the level word sits inside a noun phrase where the paper refers to its own
# genome ("This half haplotype-resolved hexaploid genome", "our tetraploid assembly").
_ANCHOR_WINDOW = 150
_SELF_REF_BEFORE_RE = re.compile(r"\b(?:this|our|the\s+present)\b[^.;]{0,60}$", re.IGNORECASE)
_SELF_REF_AFTER_RE = re.compile(r"^[^.;]{0,40}?\b(?:genome|assembly)\b", re.IGNORECASE)


def count_haploid_mentions(text: str) -> int:
    """Count organism-level uses of "haploid" in `text` (see _HAPLOID_RE)."""
    count = 0
    for match in _HAPLOID_RE.finditer(text):
        if _HAPLOID_SKIP_BEFORE_RE.search(text[max(0, match.start() - 10) : match.start()]):
            continue
        if _HAPLOID_TECHNICAL_RE.search(text[match.end() : match.end() + 30]):
            continue
        count += 1
    return count


def _is_anchored(text: str, start: int, end: int, names: list[str]) -> bool:
    """True if a target name occurs within _ANCHOR_WINDOW characters of text[start:end]
    (always True when there are no names to check against)."""
    if not names:
        return True
    if _SELF_REF_BEFORE_RE.search(text[max(0, start - 60) : start]) and _SELF_REF_AFTER_RE.search(
        text[end : end + 40]
    ):
        return True
    near = text[max(0, start - _ANCHOR_WINDOW) : end + _ANCHOR_WINDOW].lower()
    return any(re.search(rf"\b{re.escape(name)}", near) for name in names)


# Evidence weighting
_COMPOUND_BOOST = 6.0
_FORMULA_BOOST = 10.0

# Confidence cap by how many independent mentions back the winning level. A level
# named once in passing (no chromosome formula) is too weak to report confidently,
# however lopsided the vote share looks.
_SUPPORT_CONF_CAP = {1: 0.5, 2: 0.7}


def formula_matches_chromosomes(level: int, two_n: int, chromosome_count: int | None) -> bool:
    """Check a `2n = Nx = M` formula against the assembly's chromosome count.

    The assembly reports either the haploid set (M / 2), the full set (M), or a
    multiple of the base number x = M / N (subgenome / haplotype-resolved
    assemblies). Up to two extra organelle "chromosomes" (MT / Pt) are tolerated.
    A formula that fits none of these describes another organism, e.g. a
    `2n = 4x = 40` Dioscorea formula cited in a paper about A. thaliana (5)."""
    if not chromosome_count:
        return True
    base = two_n / level
    for count in (chromosome_count, chromosome_count - 1, chromosome_count - 2):
        if count <= 0:
            continue
        if count in (two_n, two_n / 2):
            return True
        if base.is_integer() and count % base == 0 and count <= two_n:
            return True
    return False


def find_level_mentions(
    text: str, window: int = 50, chromosome_count: int | None = None, names: list[str] | None = None
) -> list:
    """Find every ploidy-LEVEL mention and classify it.

    Shared by resolve_ploidy_fields (extract.py) and infer_ploidy (search.py).
    Chromosome formulas inconsistent with `chromosome_count` are dropped. A mention
    is `is_off_target` when it is a technical use of the word, or when its context
    names another species but not the target organism (`names`, see target_names).
    `is_anchored` records whether a target name occurs near the mention.
    Each item: {level, prefix, is_compound, is_ancestor, is_off_target, is_anchored,
    is_formula, context}
    """
    mentions = []

    # 2n = Nx = M  (definitive)
    for match in _CHROM_FORMULA_RE.finditer(text):
        level = int(match.group(1))
        if not (1 <= level <= 20):
            continue
        two_n = int(match.group(2))
        if not formula_matches_chromosomes(level, two_n, chromosome_count):
            logger.info(
                f"  [extract] dropped formula '{match.group(0)}': inconsistent with "
                f"{chromosome_count} assembled chromosomes"
            )
            continue
        context = text[max(0, match.start() - window) : match.end() + window]
        mentions.append(
            {
                "level": level,
                "prefix": None,
                "is_compound": True,
                "is_ancestor": False,
                "is_off_target": _names_other_species(context, names or []),
                "is_anchored": _is_anchored(text, match.start(), match.end(), names or []),
                "is_formula": True,
                "context": context,
            }
        )

    # ploidy words
    for match in _LEVEL_RE.finditer(text):
        prefix = (match.group(1) or "").lower()
        base = match.group(2).lower()
        level = _WORD_LEVEL[base]
        context = text[max(0, match.start() - window) : match.end() + window]
        is_compound = prefix in ("allo", "auto")
        is_ancestor = (not is_compound) and bool(_ANCESTOR_CONTEXT_RE.search(context))
        is_off_target = bool(_TECHNICAL_CONTEXT_RE.search(text[match.end() : match.end() + 40])) or (
            _names_other_species(context, names or [])
        )
        mentions.append(
            {
                "level": level,
                "prefix": prefix,
                "is_compound": is_compound,
                "is_ancestor": is_ancestor,
                "is_off_target": is_off_target,
                "is_anchored": _is_anchored(text, match.start(), match.end(), names or []),
                "is_formula": False,
                "context": context,
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
        match = rx.search(text)
        if not match:
            continue
        context = text[max(0, match.start() - window) : match.end() + window]
        if _MECH_ANCESTOR_RE.search(context):
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


def resolve_ploidy_fields(
    segments: list, chromosome_count: int | None = None, names: list[str] | None = None
) -> dict:
    """Resolve the 4-element ploidy tuple from weighted text segments.

    segments: iterable of (text, weight, is_priority).
    chromosome_count: NCBI assembled chromosome count, used to drop chromosome
    formulas that belong to another organism (None skips the check).
    names: target organism names (target_names); mentions that name only another
    species nearby are not counted.
    Returns the schema dict (level / mechanism / irregular / derivation + meta).
    Used by both extract.py (sections) and search.py (retrieved chunks).
    """
    level_counts: dict[int, float] = {}
    level_support: dict[int, int] = {}
    anchored_levels: set[int] = set()
    haploid_mentions = 0
    priority_derivation = None
    best_ev: dict[int, dict] = {}
    mech_votes: dict[str, float] = {}
    irregular: str | bool = False
    derivation = None
    evidence_count = 0

    for text, weight, is_priority in segments:
        if not text or weight == 0.0:
            continue

        for mention in find_level_mentions(text, chromosome_count=chromosome_count, names=names):
            if mention["is_ancestor"] or mention["is_off_target"]:
                continue

            evidence_count += 1
            # an explicit allo-/auto- term ("allotetraploid peanut") counts as two mentions
            support = 2 if mention["is_compound"] else 1
            level_support[mention["level"]] = level_support.get(mention["level"], 0) + support
            if mention["is_anchored"]:
                anchored_levels.add(mention["level"])
            if mention["is_formula"]:
                boost_weight = weight * _FORMULA_BOOST
            elif mention["is_compound"]:
                boost_weight = weight * _COMPOUND_BOOST
            else:
                boost_weight = weight
            lvl = mention["level"]
            level_counts[lvl] = level_counts.get(lvl, 0.0) + boost_weight

            rank_new = (mention["is_formula"], mention["is_compound"], is_priority)
            current = best_ev.get(lvl)
            if current is None or rank_new > current["_rank"]:
                best_ev[lvl] = {
                    "context": mention["context"],
                    "is_priority": is_priority,
                    "is_formula": mention["is_formula"],
                    "_rank": rank_new,
                }

        mechanism = detect_mechanism(text)
        if mechanism:
            mech_votes[mechanism] = mech_votes.get(mechanism, 0.0) + weight

        if not irregular:
            irregular = detect_irregular(text) or False

        if derivation is None:
            derivation = detect_derivation(text)
        if is_priority and priority_derivation is None:
            priority_derivation = detect_derivation(text)
        haploid_mentions += count_haploid_mentions(text)

    # ----- choose level -----
    if level_counts:
        level = max(level_counts, key=lambda k: level_counts[k])

        total = sum(level_counts.values())
        base_conf = level_counts[level] / total
        ev = best_ev.get(level, {})
        priority_bonus = 0.1 if ev.get("is_priority") else 0.0
        confidence = min(round(base_conf + priority_bonus, 2), 0.95)
        if not ev.get("is_formula") and level_support[level] in _SUPPORT_CONF_CAP:
            confidence = min(confidence, _SUPPORT_CONF_CAP[level_support[level]])

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
        "evidence_count": evidence_count,
        "weighted_counts": {k: round(v, 2) for k, v in level_counts.items()},
        # signals for the pipeline's precision guard
        "anchored": level is not None and level in anchored_levels,
        "haploid_mentions": haploid_mentions,
        "priority_derivation": priority_derivation,
    }


def extract_ploidy(
    weighted_sections: dict, chromosome_count: int | None = None, names: list[str] | None = None
) -> dict:
    """Extract the INSDC ploidy tuple from weighted sections by delegating to
    resolve_ploidy_fields (marking abstract/introduction/background as priority)."""
    segments = [
        (text, weight, any(p in sec.lower() for p in PRIORITY_SECTIONS))
        for sec, (text, weight) in weighted_sections.items()
    ]
    return resolve_ploidy_fields(segments, chromosome_count=chromosome_count, names=names)


def ncbi_chromosome_count(chromosome_number: dict | None) -> int | None:
    """Return the NCBI-reported chromosome count as an int, or None when it is
    missing or was only extracted from text (too unreliable to validate against)."""
    if not chromosome_number or chromosome_number.get("source") != "ncbi":
        return None
    try:
        return int(chromosome_number.get("value"))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


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

# Characters of surrounding text kept as evidence around a chromosome-number match.
CONTEXT_WINDOW = 50

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
                context_start = max(0, match.start() - CONTEXT_WINDOW)
                context_end = min(len(text), match.end() + CONTEXT_WINDOW)
                return {
                    "value": number,
                    "evidence": text[context_start:context_end],
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


def _get_ncbi_metadata(assembly: dict) -> dict:
    """Pull the NCBI-sourced assembly fields (accession, name, taxon, species,
    chromosome number and its source) that downstream extraction relies on."""
    return {
        "accession": assembly.get("assembly_accession"),
        "name": assembly.get("assembly_name"),
        "taxon_id": assembly.get("taxon_id"),
        "scientific_name": assembly.get("scientific_name"),
        "chromosome_number": assembly.get("chromosome_number"),
        "chromosome_source": assembly.get("chromosome_source"),
    }


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
    ncbi = _get_ncbi_metadata(assembly)
    ncbi_accession = ncbi["accession"]
    ncbi_name = ncbi["name"]
    ncbi_taxon = ncbi["taxon_id"]
    ncbi_species = ncbi["scientific_name"]
    ncbi_chrom = ncbi["chromosome_number"]
    ncbi_chrom_source = ncbi["chromosome_source"]

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
    names = target_names(ncbi_species, assembly.get("common_name"))
    ploidy = extract_ploidy(
        weighted_sections, chromosome_count=ncbi_chromosome_count(chromosome_number), names=names
    )
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
        "target_names": names,
        "ploidy": ploidy,
        "chromosome_number": chromosome_number,
        "cultivars": cultivars,
        "sex": sex,
    }
