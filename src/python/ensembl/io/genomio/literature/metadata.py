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
"""Assembly metadata retrieval and GoaT / NCBI-Taxonomy enrichment."""

import logging
import re as _re
from xml.etree import ElementTree as _ET

import requests

from ensembl.io.genomio.literature.ncbi import entrez_get

logger = logging.getLogger(__name__)


def _fetch_assembly_report(accession: str) -> list:
    """Query NCBI Datasets for one accession's dataset report; return the raw
    reports list (empty on any failure). Isolated so callers can retry with a
    different accession form (e.g. version-less) without duplicating HTTP logic."""
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession" f"/{accession}/dataset_report"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json().get("reports", [])
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"  [NCBI] Failed to fetch {accession}: {e}")
        return []


def fetch_assembly_metadata(accession: str) -> dict:
    """Fetch an assembly's NCBI Datasets metadata (species, taxon_id, common name,
    linked PMIDs, chromosome number), retrying without a stale version suffix.
    Returns a dict with empty/None fields if no report is found."""
    reports = _fetch_assembly_report(accession)

    # A stale version suffix (e.g. GCA_000001405.15 when .29 is current) makes
    # NCBI return an empty report. A version-less accession resolves to the
    # latest version, so retry without the ".N" suffix before giving up — this
    # keeps the pipeline working when accession lists drift out of date.
    resolved_accession = accession
    if not reports and "." in accession:
        base = accession.split(".")[0]
        logger.info(f"  [NCBI] No report for {accession}; retrying latest version ({base}) ...")
        reports = _fetch_assembly_report(base)

    if not reports:
        logger.info(f"  [NCBI] No report found for {accession}")
        return {
            "assembly_accession": accession,
            "chromosome_number": None,
            "chromosome_source": None,
        }

    report = reports[0]
    resolved_accession = report.get("accession") or resolved_accession
    assembly_info = report.get("assembly_info", {})
    organism = report.get("organism", {})
    assembly_stats = report.get("assembly_stats", {})

    chromosome_number = assembly_stats.get("total_number_of_chromosomes") or assembly_info.get(
        "chromosome_count"
    )

    raw_pmid = assembly_info.get("linked_pmid") or assembly_info.get("biosample", {}).get("linked_pmid")
    linked_pmids = raw_pmid if isinstance(raw_pmid, list) else [raw_pmid] if raw_pmid else []

    return {
        "assembly_accession": resolved_accession,
        "assembly_name": assembly_info.get("assembly_name"),
        "taxon_id": str(organism.get("tax_id") or ""),
        "scientific_name": organism.get("organism_name"),
        "common_name": organism.get("common_name"),  # NEW — mentor: search by common name too
        "linked_pmids": linked_pmids,
        "chromosome_number": chromosome_number,
        "chromosome_source": "ncbi" if chromosome_number else None,
    }


def fetch_taxonomy_common_name(taxon_id: str) -> str | None:
    """Look up a taxon's common name from NCBI Taxonomy (GenBank common name
    preferred), or None if unavailable."""
    # Caller (enrich_assembly_metadata) guards on taxon_id before calling this.
    try:
        r = entrez_get("efetch.fcgi", {"db": "taxonomy", "id": str(taxon_id)})
        root = _ET.fromstring(r.text)
        taxon = root.find("Taxon")
        if taxon is None:
            return None
        # GenbankCommonName preferred, else the first CommonName
        genbank_cname = taxon.find(".//OtherNames/GenbankCommonName")
        if genbank_cname is not None and genbank_cname.text:
            return genbank_cname.text.strip()
        common_name = taxon.find(".//OtherNames/CommonName")
        if common_name is not None and common_name.text:
            return common_name.text.strip()
    except (requests.RequestException, _ET.ParseError):
        pass
    return None


def fetch_goat_reference(taxon_id: str) -> dict:
    """Look up taxon-level reference ploidy / chromosome number from GoaT.

    GoaT stores attributes under records[0].record.attributes. Each attribute
    carries an `aggregation_source`: "direct" = a directly measured/curated
    value; "ancestor"/"descendant" = phylogenetically *inferred*. We only trust
    DIRECT values — inferred ones repeat the polyploid-ancestor trap (e.g. for
    Arabidopsis, GoaT's ploidy_inferred is 4 via the tribe, which is wrong).
    Any unexpected shape -> returns Nones so the pipeline falls back safely.
    """
    out: dict[str, int | None] = {"reference_ploidy": None, "reference_chromosome": None}
    # Caller (enrich_assembly_metadata) guards on taxon_id before calling this.
    url = "https://goat.genomehubs.org/api/v2/record"
    try:
        r = requests.get(
            url,
            params={
                "recordId": str(taxon_id),
                "result": "taxon",
                "taxonomy": "ncbi",
            },
            timeout=30,
        )
        r.raise_for_status()
        recs = r.json().get("records", [])
        if not recs:
            return out
        attrs = (recs[0].get("record", {}) or {}).get("attributes", {}) or {}

        def _direct_int(*keys: str) -> int | None:
            # try each attribute key in order; accept only directly-measured values
            for key in keys:
                attr = attrs.get(key)
                if not isinstance(attr, dict):
                    continue
                if attr.get("aggregation_source") != "direct":
                    continue  # skip inferred (ancestor/descendant) values
                for value_key in ("value", "median", "mode", "max", "min"):
                    value = attr.get(value_key)
                    if value is None:
                        continue
                    match = _re.search(r"\d+", str(value))
                    if match:
                        return int(match.group())
            return None

        out["reference_ploidy"] = _direct_int("ploidy", "ploidy_inferred")
        out["reference_chromosome"] = _direct_int("chromosome_number", "haploid_chromosome_count")
    except (requests.RequestException, ValueError, KeyError, TypeError):
        pass
    return out


def enrich_assembly_metadata(assembly: dict) -> dict:
    """Fill common_name (if missing) and attach GoaT reference values, using taxon_id."""
    taxon_id = assembly.get("taxon_id")
    if not taxon_id:
        return assembly

    if not assembly.get("common_name"):
        common_name = fetch_taxonomy_common_name(taxon_id)
        if common_name:
            assembly["common_name"] = common_name
            assembly["common_name_source"] = "ncbi_taxonomy"
            logger.info(f"      common_name filled from NCBI Taxonomy: {common_name}")

    goat = fetch_goat_reference(taxon_id)
    assembly["reference_ploidy"] = goat.get("reference_ploidy")
    assembly["reference_chromosome"] = goat.get("reference_chromosome")
    if goat.get("reference_ploidy") is not None:
        logger.info(f"      GoaT reference ploidy (taxon {taxon_id}): {goat['reference_ploidy']}")
    return assembly
