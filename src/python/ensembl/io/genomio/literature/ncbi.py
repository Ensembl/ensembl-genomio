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

import requests

from ensembl.io.genomio.literature.europepmc import search_by_pmid

logger = logging.getLogger(__name__)


# NCBI E-utilities base URL, shared by the Entrez helpers below.
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def entrez_get(endpoint: str, params: dict) -> requests.Response:
    """Call an NCBI E-utilities endpoint (e.g. "esearch.fcgi") and return the raw
    response. Raises on HTTP error; callers wrap it in their own try/except."""
    response = requests.get(f"{NCBI_EUTILS}/{endpoint}", params=params, timeout=30)  # type: ignore[arg-type]
    response.raise_for_status()
    return response


def entrez_json(endpoint: str, params: dict) -> dict:
    """Call an NCBI E-utilities endpoint and return the parsed JSON response."""
    return entrez_get(endpoint, params).json()


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


def fetch_linked_pubmed_for_assembly(accession: str, max_results: int = 5) -> list:
    """Follow the NCBI assembly -> PubMed elink to retrieve linked publications.
    Returns Europe PMC paper records ([] if no assembly UID or no links found)."""
    # Step 1: convert GCA accession to NCBI assembly UID
    try:
        data = entrez_json("esearch.fcgi", {"db": "assembly", "term": accession, "retmode": "json"})
        assembly_ids = data.get("esearchresult", {}).get("idlist", [])
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"  [NCBI Entrez elink] Assembly search failed: {e}")
        return []

    if not assembly_ids:
        logger.info(f"  [NCBI Entrez elink] No assembly ID found for {accession}")
        return []

    assembly_id = assembly_ids[0]
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
    """Mentor item 3: follow assembly -> BioProject -> PubMed to retrieve the
    canonical *reference* paper for the assembly. These are attached even when
    they carry no ploidy information (a BioProject paper still anchors species /
    common name / provenance). Pure-Entrez chain; defensive on every step."""
    # Step 1: accession -> assembly UID
    try:
        data = entrez_json("esearch.fcgi", {"db": "assembly", "term": accession, "retmode": "json"})
        asm_ids = data.get("esearchresult", {}).get("idlist", [])
    except (requests.RequestException, ValueError):
        asm_ids = []
    if not asm_ids:
        return []

    # Step 2: assembly UID -> BioProject UID(s)
    bioproject_ids = _entrez_elink("assembly", "bioproject", asm_ids[0])
    if not bioproject_ids:
        logger.info("  [BioProject] No linked BioProject found")
        return []
    logger.info(f"  [BioProject] Linked BioProject UID(s): {bioproject_ids[:5]}")

    # Step 3: BioProject UID(s) -> reference PMIDs
    pmids = []
    for bp in bioproject_ids[:3]:
        for pmid in _entrez_elink("bioproject", "pubmed", bp):
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
