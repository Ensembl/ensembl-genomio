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
"""Fetch all the sequence regions from a core database and print them in JSON format."""

__all__ = [
    "fetch_coord_systems",
    "get_seq_regions",
    "add_attribs",
    "get_synonyms",
    "get_karyotype",
]

import json
import logging
from pathlib import Path
from typing import Any, Iterator

from sqlalchemy import select
from sqlalchemy.orm import Session, joinedload

from ensembl.core.models import CoordSystem, SeqRegion, SeqRegionSynonym, SeqRegionAttrib
import ensembl.io.genomio
from ensembl.io.genomio.database import DBConnectionLite
from ensembl.io.genomio.external_db.db_map import get_external_db_map, DEFAULT_EXTERNAL_DB_MAP
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


_KARYOTYPE_STRUCTURE = {"TEL": "telomere", "ACEN": "centromere"}


def fetch_coord_systems(session: Session) -> Iterator[CoordSystem]:
    """Retrieve the coord_system metadata from the current core.

    Args:
        session: Session for the current core database.

    Yields:
        All default coord_systems in the core database.
    """
    coord_system_select = select(CoordSystem).filter(CoordSystem.attrib.like(r"%default_version%"))
    for row in session.execute(coord_system_select).unique().all():
        coord: CoordSystem = row[0]
        yield coord


def fetch_seq_regions(session: Session, coord_system: CoordSystem) -> Iterator[SeqRegion]:
    """Retrieve the seq_region metadata for a given coord_system, with accessory tables cached.

    Args:
        session: Session for the current core database.
        coord_system: Coord_system to get seq regions.

    Yields:
        All seq_regions for the coord_system.
    """
    seq_region_select = (
        select(SeqRegion)
        .where(SeqRegion.coord_system_id == coord_system.coord_system_id)
        .options(
            joinedload(SeqRegion.seq_region_synonym).joinedload(SeqRegionSynonym.external_db),
            joinedload(SeqRegion.seq_region_attrib).joinedload(SeqRegionAttrib.attrib_type),
            joinedload(SeqRegion.karyotype),
        )
    )
    for row in session.execute(seq_region_select).unique().all():
        seq_region: SeqRegion = row[0]
        yield seq_region


def get_attribs_dict(seq_region: SeqRegion) -> dict[str, Any]:
    """Returns a dict of attrib code-value for all the attributes of the given sequence region."""
    return {attrib.attrib_type.code: attrib.value for attrib in seq_region.seq_region_attrib}


def add_attribs(seq_region: dict, seq_region_attrib: dict) -> None:
    """Map seq_regions attribs to a specific name and type defined below.

    Args:
        seq_region: A seq_region dict to modify.
        seq_region_attrib: The attribs for this seq_region.
    """
    bool_attribs = {
        "circular_seq": "circular",
        "non_ref": "non_ref",
    }
    int_attribs = {
        "codon_table": "codon_table",
    }
    string_attribs = {
        "BRC4_seq_region_name": "BRC4_seq_region_name",
        "EBI_seq_region_name": "EBI_seq_region_name",
        "coord_system_tag": "coord_system_level",
        "sequence_location": "location",
    }

    for name, key in bool_attribs.items():
        # Make sure "0" means False, i.e. not added
        value = int(seq_region_attrib.get(name, "0"))
        if value:
            seq_region[key] = bool(value)

    for name, key in int_attribs.items():
        value = seq_region_attrib.get(name, "")
        if value:
            seq_region[key] = int(value)

    for name, key in string_attribs.items():
        value = seq_region_attrib.get(name, "")
        if value:
            seq_region[key] = str(value)


def get_synonyms(seq_region: SeqRegion, external_db_map: dict[str, str]) -> list[dict[str, str]]:
    """Get all synonyms for a given seq_region. Use the mapping for synonym source names.

    Args:
        seq_region: Seq_region from which the synonyms are extracted.
        external_db_map: To map the synonym source names.

    Returns:
        List of all synonyms as a dict with 'name' and 'source' keys.
    """
    synonyms = seq_region.seq_region_synonym
    syns = []
    if synonyms:
        for syn in synonyms:
            if syn.external_db:
                source = syn.external_db.db_name
                if source in external_db_map:
                    source = external_db_map[source]
                syn_obj = {"name": syn.synonym, "source": source}
            else:
                syn_obj = {"name": syn.synonym}
            syns.append(syn_obj)

    syns = sorted(syns, key=lambda syn: (syn["name"], syn.get("source", "")))
    return syns


def get_karyotype(seq_region: SeqRegion) -> list[dict[str, str]]:
    """Given a seq_region, extract the karyotype bands.

    Args:
        seq_region: The seq_region from which the karyotype bands are extracted.

    Returns:
        List of all karyotype bands as a dict with values 'start', 'end', 'name' 'stain', 'structure'.
    """
    bands = seq_region.karyotype
    kars = []
    if bands:
        for band in bands:
            kar = {"start": band.seq_region_start, "end": band.seq_region_end}
            if band.band:
                kar["name"] = band.band
            if band.stain:
                kar["stain"] = band.stain
                structure = _KARYOTYPE_STRUCTURE.get(band.stain, "")
                if structure:
                    kar["structure"] = structure
            kars.append(kar)

    kars = sorted(kars, key=lambda kar: kar.get("name", ""))
    return kars


def get_added_sequence(seq_region: SeqRegion) -> dict[str, str | dict[str, str]]:
    """Extracts added sequence information of the given sequence region.

    Args:
        seq_region: Sequence region.

    Returns:
        Accession as well as assembly and annotation provider information of the added sequence.
    """
    attribs = get_attribs_dict(seq_region)
    accession = attribs.get("added_seq_accession")
    if accession is None:
        return {}

    added_sequence = {
        "accession": accession,
    }

    # Assembly provider
    assembly_provider = {
        "name": attribs.get("added_seq_asm_pr_nam"),
        "url": attribs.get("added_seq_asm_pr_url"),
    }
    if assembly_provider["name"] and assembly_provider["url"]:
        added_sequence["assembly_provider"] = assembly_provider

    # annotation provider
    annotation_provider = {
        "name": attribs.get("added_seq_ann_pr_nam"),
        "url": attribs.get("added_seq_ann_pr_url"),
    }
    if annotation_provider["name"] and annotation_provider["url"]:
        added_sequence["annotation_provider"] = annotation_provider

    return added_sequence


def get_seq_regions(session: Session, external_db_map: dict) -> list[SeqRegion]:
    """Returns all the sequence regions from the current core database.

    Include synonyms, attribs and karyotypes. Only the top level sequences are exported.

    Args:
        session: Session from the current core database.
        external_db_map: Mapping of external_db names for the synonyms.

    """
    seq_regions = []

    for coord_system in fetch_coord_systems(session):
        logging.debug(f"Dump coord {coord_system.name}")
        for seqr in fetch_seq_regions(session, coord_system):
            seq_region: dict[str, Any] = {}
            seq_region = {"name": seqr.name, "length": seqr.length}
            synonyms = get_synonyms(seqr, external_db_map)
            if synonyms:
                seq_region["synonyms"] = synonyms

            attribs = get_attribs_dict(seqr)
            if not attribs or "toplevel" not in attribs:
                # Skip seq_region without attribs or not toplevel
                continue
            add_attribs(seq_region, attribs)

            karyotype = get_karyotype(seqr)
            if karyotype:
                seq_region["karyotype_bands"] = karyotype

            added_seq = get_added_sequence(seqr)
            if added_seq:
                seq_region["added_sequence"] = added_seq

            if "coord_system_level" not in seq_region:
                seq_region["coord_system_level"] = coord_system.name

            seq_regions.append(seq_region)

    seq_regions = sorted(seq_regions, key=lambda seqr: (seqr["coord_system_level"], seqr["name"]))
    return seq_regions


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Fetch all the sequence regions from a core database and print them in JSON format."
    )
    parser.add_server_arguments(include_database=True)
    parser.add_argument_src_path(
        "--external_db_map", default=DEFAULT_EXTERNAL_DB_MAP.resolve(), help="File with external_db mapping"
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnectionLite(args.url)

    external_map_path = Path(args.external_db_map)
    external_map = get_external_db_map(external_map_path)

    with dbc.session_scope() as session:
        seq_regions = get_seq_regions(session, external_map)

    print(json.dumps(seq_regions, indent=2, sort_keys=True))
