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
    "MapFormatError",
    "get_external_db_map",
    "get_coord_systems",
    "get_seq_regions",
    "add_attribs",
    "get_synonyms",
    "get_karyotype",
]

import json
import logging
from pathlib import Path
from typing import Any, Dict, List

from sqlalchemy import select
from sqlalchemy.orm import Session, joinedload

from ensembl.core.models import CoordSystem, SeqRegion, SeqRegionSynonym, SeqRegionAttrib
from ensembl.io.genomio.database import DBConnectionLite
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


_ROOT_DIR = Path(__file__).parent / "../../../../../.."
_DEFAULT_MAP = _ROOT_DIR / "config/external_db_map/default.txt"
_KARYOTYPE_STRUCTURE = {"TEL": "telomere", "ACEN": "centromere"}


class MapFormatError(Exception):
    """Error when parsing the db map file."""


def get_external_db_map(map_file: Path) -> Dict:
    """Class method, set up the map for all SeqRegion objects"""
    db_map: Dict[str, str] = {}
    with map_file.open("r") as map_fh:
        for line in map_fh:
            line = line.rstrip()
            if line.startswith("#") or line.startswith(" ") or line == "":
                continue
            parts = line.split("\t")
            if not parts[0] or not parts[1]:
                raise MapFormatError(f"External db file is not formatted correctly for: {line}")
            db_map[parts[1]] = parts[0]
    return db_map


def get_coord_systems(session: Session) -> List[CoordSystem]:
    """Retrieve the coord_system metadata from the current core.

    Args:
        session (Session): Session for the current core.

    Returns:
        List[CoordSystem]: All coord_systems in the core.
    """
    coord_systems: List[CoordSystem] = []
    coord_stmt = select(CoordSystem).filter(CoordSystem.attrib.like("%default_version%"))
    for row in session.execute(coord_stmt).unique().all():
        coord_systems.append(row[0])
    return coord_systems


def get_seq_regions(session: Session, external_db_map: dict) -> List[SeqRegion]:
    """Returns all the sequence regions from the current core database.

    Include synonyms, attribs and karyotypes. Only the top level sequences are exported.

    Args:
        session (Session): Session from the current core.
        external_db_map (dict): Mapping of external_db names for the synonyms.

    """
    coord_systems = get_coord_systems(session)
    seq_regions = []

    for coord_system in coord_systems:
        logging.debug(f"Dump coord {coord_system.name}")
        seqr_stmt = (
            select(SeqRegion)
            .where(SeqRegion.coord_system_id == coord_system.coord_system_id)
            .options(
                joinedload(SeqRegion.seq_region_synonym).joinedload(SeqRegionSynonym.external_db),
                joinedload(SeqRegion.seq_region_attrib).joinedload(SeqRegionAttrib.attrib_type),
                joinedload(SeqRegion.karyotype),
            )
        )
        for row in session.execute(seqr_stmt).unique().all():
            seqr: SeqRegion = row[0]
            seq_region: Dict[str, Any] = {}
            seq_region = {"name": seqr.name, "length": seqr.length}
            synonyms = get_synonyms(seqr, external_db_map)
            if synonyms:
                seq_region["synonyms"] = synonyms

            attribs = get_attribs_dict(seqr)
            if attribs:
                if "toplevel" not in attribs:
                    continue
                add_attribs(seq_region, attribs)
            else:
                # Skip seq_region without attribs, not toplevel
                continue

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


def add_attribs(seq_region: Dict, attrib_dict: Dict) -> None:
    """Map seq_regions attribs to a specific name and type defined below.

    Args:
        seq_region (Dict): A seq_region dict to modify.
        attrib_dict (Dict): The attribs for this seq_region.
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
        value = attrib_dict.get(name)
        if value:
            seq_region[key] = bool(value)

    for name, key in int_attribs.items():
        value = attrib_dict.get(name)
        if value:
            seq_region[key] = int(value)

    for name, key in string_attribs.items():
        value = attrib_dict.get(name)
        if value:
            seq_region[key] = str(value)


def get_synonyms(seq_region: SeqRegion, external_db_map: dict) -> List:
    """Get all synonyms for a given seq_region. Use the mapping for synonym source names.

    Args:
        seq_region (SeqRegion): Seq_region from which the synonyms are extracted.
        external_db_map (dict): To map the synonym source names.

    Returns:
        List: All synonyms as a dict with 'name' and 'source' keys.
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


def _get_attribs(seq_region: SeqRegion) -> List:
    """Given a seq_region, extract the attribs as value-source items.

    Args:
        seq_region (SeqRegion): The seq_region from which the attribs are extracted.

    Returns:
        All attributes as a list of dictionaries with 'value' and 'source' keys.
    """
    attribs = seq_region.seq_region_attrib
    atts = []
    if attribs:
        for attrib in attribs:
            att_obj = {"value": attrib.value, "source": attrib.attrib_type.code}
            atts.append(att_obj)
    return atts


def get_attribs_dict(seq_region: SeqRegion) -> Dict[str, Any]:
    """Extracts all the attributes of the given sequence region.

    Args:
        seq_region: Sequence region.

    Returns:
        Pairs of source and value for each attribute.
    """

    attribs = _get_attribs(seq_region)
    attrib_dict = {attrib["source"]: attrib["value"] for attrib in attribs}
    return attrib_dict


def get_karyotype(seq_region: SeqRegion) -> List:
    """Given a seq_region, extract the karyotype bands.

    Args:
        seq_region (SeqRegion): The seq_region from which the karyotype bands are extracted.

    Returns:
        List: All karyotype bands as a dict with values 'start', 'end', 'name' 'stain', 'structure'.
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

    kars = sorted(kars, key=lambda kar: kar["name"])
    return kars


def get_added_sequence(seq_region: SeqRegion) -> Dict:
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


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Fetch all the sequence regions from a core database and print them in JSON format."
    )
    parser.add_server_arguments(include_database=True)
    parser.add_argument_src_path(
        "--external_db_map", default=_DEFAULT_MAP.resolve(), help="File with external_db mapping"
    )
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnectionLite(args.url)

    external_map_path = Path(args.external_db_map)
    external_map = get_external_db_map(external_map_path)

    with dbc.session_scope() as session:
        seq_regions = get_seq_regions(session, external_map)

    print(json.dumps(seq_regions, indent=2, sort_keys=True))
