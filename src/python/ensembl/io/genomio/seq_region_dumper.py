#!/usr/bin/env python
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
"""Generates one JSON file per metadata type inside `manifest`, including the manifest itself.

Can be imported as a module and called as a script as well, with the same parameters and expected outcome.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import argschema
from sqlalchemy import select
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session, joinedload

from ensembl.database import DBConnection
from ensembl.core.models import CoordSystem, SeqRegion, SeqRegionSynonym, SeqRegionAttrib

ROOT_DIR = Path(__file__).parent / "../../../../.."
DEFAULT_MAP = ROOT_DIR / "config/external_db_map/default.txt"
KARYOTYPE_STRUCTURE = {"TEL": "telomere", "ACEN": "centromere"}


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
    """Retrieve the seq_region metadata from the current core.
    Include synonyms, attribs and karyotypes.
    Only the top level sequences are exported.

    Args:
        session (Session): Session from the current core.
        external_db_map (dict): Mapping of external_db names for the synonyms.

    Returns:
        List[SeqRegion]: All seq_regions in the core.
    """
    coord_systems = get_coord_systems(session)
    seq_regions = []

    for coord_system in coord_systems:
        print(f"Dump coord {coord_system.name}")
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

            attribs = get_attribs(seqr)
            if attribs:
                attrib_dict = {attrib["source"]: attrib["value"] for attrib in attribs}
                if "toplevel" not in attrib_dict:
                    continue
                add_attribs(seq_region, attrib_dict)
            else:
                # Skip seq_region without attribs, not toplevel
                continue

            karyotype = get_karyotype(seqr)
            if karyotype:
                seq_region["karyotype_bands"] = karyotype

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


def get_attribs(seq_region: SeqRegion) -> List:
    """Given a seq_region, extract the attribs as value-source items.

    Args:
        seq_region (SeqRegion): The seq_region from which the attribs are extracted.

    Returns:
        List: All attribs as a dict with 'value' and 'source' keys.
    """
    attribs = seq_region.seq_region_attrib
    atts = []
    if attribs:
        for attrib in attribs:
            att_obj = {"value": attrib.value, "source": attrib.attrib_type.code}
            atts.append(att_obj)
    return atts


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
                structure = KARYOTYPE_STRUCTURE.get(band.stain, "")
                if structure:
                    kar["structure"] = structure
            kars.append(kar)

    kars = sorted(kars, key=lambda kar: kar["name"])
    return kars


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    # Server parameters
    host = argschema.fields.String(
        required=True, metadata={"description": "Host to the server with EnsEMBL databases"}
    )
    port = argschema.fields.Integer(required=True, metadata={"description": "Port to use"})
    user = argschema.fields.String(required=True, metadata={"description": "User to use"})
    password = argschema.fields.String(required=False, metadata={"description": "Password to use"})
    database = argschema.fields.String(required=True, metadata={"description": "Database to use"})
    external_db_map = argschema.fields.files.InputFile(
        required=False,
        default=str(DEFAULT_MAP),
        metadata={"description": "File with external_db mapping"},
    )


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    args = mod.args

    db_url = URL.create(
        "mysql",
        mod.args["user"],
        mod.args.get("password"),
        mod.args["host"],
        mod.args["port"],
        mod.args.get("database"),
    )
    dbc = DBConnection(db_url)

    external_map_path = Path(mod.args.get("external_db_map"))
    external_map = get_external_db_map(external_map_path)

    with dbc.session_scope() as session:
        seq_regions = get_seq_regions(session, external_map)

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(seq_regions, indent=2, sort_keys=True))
    else:
        print(seq_regions)


if __name__ == "__main__":
    main()
