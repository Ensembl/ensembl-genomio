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
"""Europe PMC search and article/supplementary text retrieval."""

import io
import logging
import re as _re
import time
import zipfile

import requests

logger = logging.getLogger(__name__)


def search_europe_pmc(query: str, max_results: int = 5, retries: int = 3) -> list:
    """Search Europe PMC for `query` and return up to `max_results` paper records,
    retrying transient request failures up to `retries` times."""
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": query,
        "format": "json",
        "pageSize": max_results,
        "resultType": "core",
    }
    for attempt in range(retries):
        try:
            response = requests.get(url, params=params, timeout=30)  # type: ignore[arg-type]
            response.raise_for_status()
            return response.json().get("resultList", {}).get("result", [])
        except requests.RequestException as e:
            if attempt < retries - 1:
                wait = 2**attempt
                logger.warning(f"  [EuropePMC] Attempt {attempt+1} failed, retrying in {wait}s ...")
                time.sleep(wait)
            else:
                logger.warning(f"  [EuropePMC] All {retries} attempts failed: {e}")
                return []
    return []


def search_by_pmid(pmid: str) -> dict | None:
    """Fetch one paper's Europe PMC record by PMID, or None if not found."""
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": f"EXT_ID:{pmid} AND SRC:MED",
        "format": "json",
        "resultType": "core",
    }
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        results = response.json().get("resultList", {}).get("result", [])
        return results[0] if results else None
    except (requests.RequestException, ValueError):
        return None


def search_europe_pmc_by_name(
    taxon_id: str,
    scientific_name: str,
    common_name: str = "",
    max_results: int = 5,
) -> list:
    """Search Europe PMC for genome papers by scientific/common name, escalating
    from full name to binomial to genus to taxon-ID queries and keeping only
    genus-relevant hits. Returns up to max_results records ([] if none)."""
    tokens = scientific_name.split()  # "".split() -> [], so no guard needed
    genus = tokens[0] if tokens else ""
    # Binomial = genus + species epithet. For a subspecies / variety / hybrid
    # (e.g. "Oryza meyeriana var. indandamanica") the full trinomial rarely
    # appears in a paper, but the binomial ("Oryza meyeriana") usually does.
    binomial = " ".join(tokens[:2])

    # Multi-identifier name clause: scientific name, plus common name when
    # it is specific enough to be useful (>= 4 chars avoids junk like "fly").
    name_clause = f'"{scientific_name}"'
    if common_name and len(common_name) >= 4:
        name_clause = f'("{scientific_name}" OR "{common_name}")'

    queries = [
        # 1st priority: name(s) + genome assembly
        f'{name_clause} AND "genome assembly"',
        # 2nd priority: name(s) + genome sequence / ploidy / chromosome
        f'{name_clause} AND ("genome sequence" OR ploidy OR chromosome)',
        # 3rd priority: scientific name + genome sequence
        f'"{scientific_name}" AND "genome sequence"',
    ]

    # 4th priority: taxonomic fallback. When the full name is a trinomial
    # (subspecies / variety / hybrid), drop to the binomial, then the genus,
    # before the taxon-ID last resort. The genus relevance filter below still
    # applies, so this widens recall without inviting unrelated papers.
    if binomial and binomial != scientific_name:
        queries += [
            f'"{binomial}" AND "genome assembly"',
            f'"{binomial}" AND ("genome sequence" OR ploidy OR chromosome)',
        ]

    if taxon_id:
        queries += [
            # 5th priority: taxon ID fallback (last resort)
            f"TAXONOMY:{taxon_id} AND genome assembly",
        ]

    for query in queries:
        results = search_europe_pmc(query, max_results=max_results * 3)
        if not results:
            continue

        # Keep only results that actually mention the genus — guards against
        # false positives from an ambiguous common name.
        relevant = [
            r for r in results if genus.lower() in (r.get("title", "") + r.get("abstractText", "")).lower()
        ]

        if relevant:
            logger.info(f"  [EuropePMC] Query: {query}")
            logger.info(
                f"              Found {len(results)} results, "
                f"{len(relevant)} relevant to {scientific_name}"
            )
            return relevant[:max_results]

    logger.info(f"  [EuropePMC] No relevant papers found for {scientific_name}")
    return []


def search_europe_pmc_fulltext(scientific_name: str, max_results: int = 5) -> list:
    """Find candidate source papers that mention the exact species binomial in
    their FULL TEXT even when the title/abstract does not — e.g. a multi-species
    or methods paper that assembled this genome (the true source of many
    assemblies). Europe PMC's search covers open-access full text, so a quoted
    binomial query surfaces them. We deliberately skip the genus-in-title guard
    here (that guard drops exactly these papers) and instead return only
    open-access hits, tagged for a full-text confirmation pass in the caller."""
    if not scientific_name:
        return []
    query = f'"{scientific_name}" AND (genome OR assembly OR chromosome OR ploidy)'
    candidates = []
    for paper in search_europe_pmc(query, max_results=max_results * 3):
        if paper.get("pmcid") or paper.get("pmCid"):
            paper["_needs_fulltext_confirm"] = True
            candidates.append(paper)
        if len(candidates) >= max_results * 2:
            break
    if candidates:
        logger.info(f"  [EuropePMC] Full-text species search: {len(candidates)} open-access candidate(s)")
    return candidates


def validate_paper_ids(paper: dict) -> dict:
    """Check a paper's PMCID / PMID against Europe PMC and report which text tier
    (full text, abstract, or none) is actually retrievable for it."""
    pmcid = paper.get("pmcid")
    pmid = paper.get("pmid")
    doi = paper.get("doi")

    pmcid_valid = False
    if pmcid:
        try:
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
            r = requests.head(url, timeout=10)
            pmcid_valid = r.status_code == 200
        except requests.RequestException:
            pmcid_valid = False

    pmid_valid = False
    if pmid:
        try:
            result = search_by_pmid(pmid)
            pmid_valid = result is not None
        except requests.RequestException:
            pmid_valid = False

    return {
        "pmcid": pmcid,
        "pmid": pmid,
        "doi": doi,
        "pmcid_valid": pmcid_valid,
        "pmid_valid": pmid_valid,
        "doi_present": bool(doi),
        "has_fulltext": pmcid_valid,
        "best_available": ("fulltext" if pmcid_valid else "abstract" if pmid_valid else "none"),
    }


def get_fulltext_by_pmcid(pmcid: str) -> str | None:
    """Download an article's full-text XML from Europe PMC by PMCID, or None."""
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.text
    except requests.RequestException:
        pass
    return None


def _strip_markup(s: str) -> str:
    try:
        from bs4 import BeautifulSoup

        return BeautifulSoup(s, "html.parser").get_text(" ")
    except Exception:
        return _re.sub(r"<[^>]+>", " ", s)


def _read_zip(raw: bytes, member: str) -> str:
    """Read one member from an in-memory zip (xlsx/docx are zip archives) and
    return its decoded text, or '' if the archive or the member is unreadable."""
    try:
        archive = zipfile.ZipFile(io.BytesIO(raw))
        if member not in archive.namelist():
            return ""
        return archive.read(member).decode("utf-8", errors="ignore")
    except Exception:
        return ""


def _xlsx_shared_strings(raw: bytes) -> str:
    """Extract shared-string text from an .xlsx (itself a zip) without openpyxl.
    This catches column headers / row labels ('2n', 'chromosome number', ploidy
    words) that carry the ploidy signal; numeric cell values are not needed for
    that. Returns '' if the file is not a readable xlsx."""
    xml = _read_zip(raw, "xl/sharedStrings.xml")
    return " ".join(_re.findall(r"<t[^>]*>([^<]+)</t>", xml))


def _docx_text(raw: bytes) -> str:
    """Extract paragraph text from a .docx (itself a zip) without python-docx:
    read word/document.xml and pull the <w:t> text runs. Supplementary tables
    are very commonly distributed as .docx, so this materially widens coverage.
    Returns '' if the file is not a readable docx."""
    xml = _read_zip(raw, "word/document.xml")
    return " ".join(_re.findall(r"<w:t[^>]*>([^<]*)</w:t>", xml))


def get_supplementary_text_by_pmcid(pmcid: str, max_chars: int = 40000) -> str | None:
    """Download the Europe PMC supplementary-files ZIP for a PMC article and
    extract plain text from the members we can read without heavy dependencies:
    txt/csv/tsv, html/xml, and xlsx shared strings. Ploidy evidence —
    chromosome counts, karyotype tables, '2n = ...' formulas — often lives in
    supplementary tables rather than the article body, so this widens what the
    extractor can see. Binary formats we cannot read dependency-free (pdf, docx,
    images) are skipped. Returns one concatenated string, or None if nothing
    usable was found."""
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
    try:
        resp = requests.get(url, timeout=60)
        if resp.status_code != 200 or not resp.content:
            return None
    except requests.RequestException:
        return None

    try:
        zf = zipfile.ZipFile(io.BytesIO(resp.content))
    except zipfile.BadZipFile:
        return None

    parts, total = [], 0
    for name in zf.namelist():
        if total >= max_chars:
            break
        lower = name.lower()
        try:
            raw = zf.read(name)
        except Exception:
            continue
        if lower.endswith((".txt", ".csv", ".tsv", ".tab")):
            text = raw.decode("utf-8", errors="ignore")
        elif lower.endswith((".html", ".htm", ".xml")):
            text = _strip_markup(raw.decode("utf-8", errors="ignore"))
        elif lower.endswith(".xlsx"):
            text = _xlsx_shared_strings(raw)
        elif lower.endswith(".docx"):
            text = _docx_text(raw)
        else:
            continue  # pdf / legacy .doc / images: not readable dependency-free
        text = _re.sub(r"[ \t\r\n\f\v]+", " ", text).strip()
        if text:
            snippet = text[: max_chars - total]
            parts.append(f"[Supplementary file: {name}]\n{snippet}")
            total += len(snippet)

    if not parts:
        return None
    logger.info(f"  [supp] extracted text from {len(parts)} supplementary file(s) for {pmcid}")
    return "\n\n".join(parts)


def get_abstract_by_pmid(pmid: str) -> str | None:
    """Return the abstract text for a PMID via Europe PMC, or None."""
    paper = search_by_pmid(pmid)
    if paper:
        return paper.get("abstractText")
    return None


def get_best_available_text(paper: dict) -> dict:
    """Return the best retrievable text for a paper as {source, text, ...}: full
    text (+ supplementary) by PMCID, else abstract by PMID, else the search-result
    abstract, else source 'none'."""
    pmcid = paper.get("pmcid")
    pmid = paper.get("pmid")

    if pmcid:
        fulltext = get_fulltext_by_pmcid(pmcid)
        if fulltext:
            out = {"source": "fulltext", "pmcid": pmcid, "text": fulltext}
            supp = get_supplementary_text_by_pmcid(pmcid)
            if supp:
                out["supplementary"] = supp
            return out

    if pmid:
        abstract = get_abstract_by_pmid(pmid)
        if abstract:
            return {"source": "abstract", "pmid": pmid, "text": abstract}

    abstract = paper.get("abstractText")
    if abstract:
        return {"source": "search_result", "text": abstract}

    return {"source": "none", "text": None}
