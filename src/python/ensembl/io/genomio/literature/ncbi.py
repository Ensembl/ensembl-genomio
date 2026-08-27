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
"""NCBI E-utilities (Entrez) helpers and NCBI-linked paper retrieval."""

import logging
import os
import time
from functools import lru_cache

import requests

from ensembl.io.genomio.literature.europepmc import search_by_pmid

logger = logging.getLogger(__name__)


# NCBI E-utilities base URL, shared by the Entrez helpers below.
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Optional NCBI API key: raises the E-utilities rate limit from 3 req/s to
# 10 req/s. Set as an env var; falls back to unauthenticated calls if unset.
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")


def entrez_get(endpoint: str, params: dict, retries: int = 3) -> requests.Response:
    """Call an NCBI E-utilities endpoint (e.g. "esearch.fcgi") and return the raw
    response, retrying transient request failures with exponential backoff and
    adding the NCBI API key when configured. Raises the last error if all
    retries fail; callers wrap it in their own try/except."""
    if NCBI_API_KEY:
        params = {**params, "api_key": NCBI_API_KEY}
    last_exc: requests.RequestException | None = None
    for attempt in range(retries):
        try:
            response = requests.get(f"{NCBI_EUTILS}/{endpoint}", params=params, timeout=30)  # type: ignore[arg-type]
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            last_exc = e
            if attempt < retries - 1:
                wait = 2**attempt
                logger.warning(
                    f"  [Entrez] {endpoint} attempt {attempt + 1}/{retries} failed ({e}); "
                    f"retrying in {wait}s ..."
                )
                time.sleep(wait)
            else:
                logger.warning(f"  [Entrez] {endpoint} failed after {retries} attempts: {e}")
    assert last_exc is not None
    raise last_exc


def entrez_json(endpoint: str, params: dict) -> dict:
    """Call an NCBI E-utilities endpoint and return the parsed JSON response."""
    return entrez_get(endpoint, params).json()


@lru_cache(maxsize=256)
def resolve_assembly_uid(accession: str) -> str | None:
    """accession -> NCBI assembly UID, cached per accession so
    fetch_linked_pubmed_for_assembly and fetch_bioproject_reference_papers
    (which both need this) don't each re-issue the same esearch call, and
    repeated lookups for the same accession within one process are free."""
    try:
        data = entrez_json("esearch.fcgi", {"db": "assembly", "term": accession, "retmode": "json"})
        ids = data.get("esearchresult", {}).get("idlist", [])
        return ids[0] if ids else None
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"  [NCBI Entrez] Assembly search failed for {accession}: {e}")
        return None


def _entrez_elink(dbfrom: str, db: str, uid: str) -> list:
    """Return NCBI elink links from `dbfrom` to `db` for one UID ([] on any failure)."""
    try:
        data = entrez_json("elink.fcgi", {"dbfrom": dbfrom, "db": db, "id": uid, "retmode": "json"})
        linksets = data.get("linksets", [{}])[0]
        for ldb in linksets.get("linksetdbs", []):
            if ldb.get("dbto") == db:
                return ldb.get("links", [])
    except (requests.RequestException, ValueError, KeyError, IndexError):
        pass
    return []


def _entrez_elink_batch(dbfrom: str, db: str, uids: list) -> dict:
    """elink for MULTIPLE uids in a single HTTP call. NCBI's elink accepts
    repeated `id` params and returns one linkset per input uid, so this
    replaces an N-call loop (one elink per uid) with one call regardless of
    how many uids are passed. Returns {uid: [linked_ids]} ({} on failure)."""
    if not uids:
        return {}
    try:
        response = entrez_get("elink.fcgi", {"dbfrom": dbfrom, "db": db, "id": uids, "retmode": "json"})
        data = response.json()
        out: dict[str, list] = {}
        for linkset in data.get("linksets", []):
            ids = linkset.get("ids") or []
            if not ids:
                continue
            uid = str(ids[0])
            for ldb in linkset.get("linksetdbs", []):
                if ldb.get("dbto") == db:
                    out[uid] = ldb.get("links", [])
        return out
    except (requests.RequestException, ValueError, KeyError, IndexError) as e:
        logger.warning(f"  [Entrez] batched elink {dbfrom}->{db} failed: {e}")
        return {}


def fetch_linked_pubmed_for_assembly(accession: str, max_results: int = 5) -> list:
    """Follow the NCBI assembly -> PubMed elink to retrieve linked publications.
    Returns Europe PMC paper records ([] if no assembly UID or no links found)."""
    # Step 1: convert GCA accession to NCBI assembly UID (cached)
    assembly_id = resolve_assembly_uid(accession)
    if not assembly_id:
        logger.info(f"  [NCBI Entrez elink] No assembly ID found for {accession}")
        return []

    logger.info(f"  [NCBI Entrez elink] Assembly ID: {assembly_id}")

    # Step 2: find linked PubMed articles via elink
    pmids = _entrez_elink("assembly", "pubmed", assembly_id)

    if not pmids:
        logger.info(f"  [NCBI Entrez elink] No linked publications found")
        return []

    logger.info(f"  [NCBI Entrez elink] Found {len(pmids)} linked PMID(s): {pmids[:5]}")

    # Step 3: fetch paper details from Europe PMC using PMIDs
    papers = []
    for pmid in pmids[:max_results]:
        paper = search_by_pmid(str(pmid))
        if paper:
            papers.append(paper)
            logger.info(f"  [NCBI Entrez elink] {pmid}: {paper.get('title', '')[:60]}...")

    return papers


def fetch_bioproject_reference_papers(accession: str, max_results: int = 5) -> list:
    """Follow assembly -> BioProject -> PubMed to retrieve the canonical *reference*
    paper for the assembly. These are attached even when they carry no ploidy
    information (a BioProject paper still anchors species / common name /
    provenance). Pure-Entrez chain; defensive on every step."""
    # Step 1: accession -> assembly UID (cached)
    assembly_id = resolve_assembly_uid(accession)
    if not assembly_id:
        return []

    # Step 2: assembly UID -> BioProject UID(s)
    bioproject_ids = _entrez_elink("assembly", "bioproject", assembly_id)
    if not bioproject_ids:
        logger.info("  [BioProject] No linked BioProject found")
        return []
    logger.info(f"  [BioProject] Linked BioProject UID(s): {bioproject_ids[:5]}")

    # Step 3: BioProject UID(s) -> reference PMIDs (one batched elink call)
    target_bps = bioproject_ids[:3]
    links_by_bp = _entrez_elink_batch("bioproject", "pubmed", target_bps)
    pmids: list = []
    for bp in target_bps:
        for pmid in links_by_bp.get(bp, []):
            if pmid not in pmids:
                pmids.append(pmid)
    if not pmids:
        logger.info("  [BioProject] No reference publication linked to BioProject")
        return []
    logger.info(f"  [BioProject] Found {len(pmids)} reference PMID(s): {pmids[:5]}")

    # Step 4: fetch paper details, tag as reference paper
    papers = []
    for pmid in pmids[:max_results]:
        paper = search_by_pmid(str(pmid))
        if paper:
            paper["retrieval_source"] = "bioproject"
            paper["is_reference_paper"] = True
            papers.append(paper)
            logger.info(f"  [BioProject] {pmid}: {paper.get('title', '')[:60]}...")
    return papers


def search_ncbi_entrez(scientific_name: str, max_results: int = 5) -> list:
    """NCBI Entrez PMC search fallback: query PubMed Central for open-access
    genome-assembly papers by scientific name. Used when accession links and the
    Europe PMC name search fail to surface a strong genome paper."""
    logger.info(f"  [NCBI Entrez] Searching PubMed Central for {scientific_name} ...")

    # Step 1: search for PMC IDs
    try:
        data = entrez_json(
            "esearch.fcgi",
            {
                "db": "pmc",
                "term": f'"{scientific_name}" AND "genome assembly" AND "open access"[filter]',
                "retmax": max_results * 2,
                "retmode": "json",
            },
        )
        pmc_ids = data.get("esearchresult", {}).get("idlist", [])
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"  [NCBI Entrez] Search failed: {e}")
        return []

    if not pmc_ids:
        logger.info(f"  [NCBI Entrez] No PMC IDs found")
        return []

    # Step 2: fetch summaries
    try:
        data = entrez_json(
            "esummary.fcgi",
            {"db": "pmc", "id": ",".join(pmc_ids[:max_results]), "retmode": "json"},
        )
        summaries = data.get("result", {})
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"  [NCBI Entrez] Summary fetch failed: {e}")
        return []

    # Step 3: convert to Europe PMC format
    papers = []
    for pmc_id in pmc_ids[:max_results]:
        summary = summaries.get(pmc_id, {})
        if not summary:
            continue
        pmcid = f"PMC{pmc_id}"
        papers.append(
            {
                "title": summary.get("title", ""),
                "pmcid": pmcid,
                "pmid": summary.get("pmid", ""),
                "doi": summary.get("doi", ""),
                "abstractText": "",
                "source": "ncbi_entrez",
            }
        )
        logger.info(f"  [NCBI Entrez] {pmcid}: {summary.get('title', '')[:60]}...")

    return papers
