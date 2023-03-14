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
from typing import Dict, List

import argschema
from sqlalchemy import select
from sqlalchemy.orm import Session

from ensembl.database import DBConnection
from ensembl.core.models import SeqRegion

ROOT_DIR = Path(__file__).parent / "../../../../.."
DEFAULT_MAP = ROOT_DIR / "data/external_db_map_default.txt"


def get_external_db_map(map_file: Path) -> Dict:
    """Class method, set up the map for all SeqRegion objects"""
    db_map = dict()
    with map_file.open("r") as map_fh:
        for line in map_fh:
            line = line.rstrip()
            if line.startswith("#") or line.startswith(" ") or line == "":
                continue
            parts = line.split("\t")
            if not parts[0] or not parts[1]:
                raise Exception(f"External db file is not formatted correctly for: {line}")
            else:
                db_map[parts[1]] = parts[0]
    return db_map


def get_seq_regions(session: Session, external_db_map: dict) -> List[SeqRegion]:
    seqr_stmt = select(SeqRegion)
    seq_regions = []
    for row in session.execute(seqr_stmt).unique().all():
        seqr: SeqRegion = row[0]
        seq_region = dict()
        seq_region = {
            "name": seqr.name,
            "length": seqr.length
        }
        synonyms = get_synonyms(seqr, external_db_map)
        if synonyms:
            seq_region["synonyms"] = synonyms

        attribs = get_attribs(seqr)
        if attribs:
            seq_region["attribs"] = attribs

        karyotype = get_karyotype(seqr)
        if karyotype:
            seq_region["karyotype"] = karyotype

        seq_regions.append(seq_region)

    return seq_regions


def get_synonyms(seq_region: SeqRegion, external_db_map: dict) -> List:
    synonyms = seq_region.seq_region_synonym
    syns = []
    if synonyms:
        for syn in synonyms:
            source = syn.external_db.db_name
            if source in external_db_map:
                source = external_db_map[source]
            syn_obj = {
                "synonym": syn.synonym,
                "source": source
            }
            syns.append(syn_obj)
    return syns


def get_attribs(seq_region: SeqRegion) -> List:
    attribs = seq_region.seq_region_attrib
    atts = []
    if attribs:
        for attrib in attribs:
            att_obj = {
                "value": attrib.value,
                "source": attrib.attrib_type.code
            }
            atts.append(att_obj)
    return atts


def get_karyotype(seq_region: SeqRegion) -> List:
    bands = seq_region.karyotype
    kars = []
    if bands:
        for band in bands:
            kar = {
                "start": band.seq_region_start,
                "end": band.seq_region_end
            }
            if band.band:
                kar["band"] = band.band
            if band.stain:
                kar["stain"] = band.stain
            kars.append(kar)
    return kars


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    # Server parameters
    host = argschema.fields.String(required=True, metadata={
        "description": "Host to the server with EnsEMBL databases"
    })
    port = argschema.fields.Integer(required=True, metadata={
        "description": "Port to use"
    })
    user = argschema.fields.String(required=True, metadata={
        "description": "User to use"
    })
    password = argschema.fields.String(required=False, metadata={
        "description": "Password to use"
    })
    database = argschema.fields.String(required=True, metadata={
        "description": "Database to use"
    })
    external_db_map = argschema.fields.files.InputFile(
        required=False,
        dump_default=str(DEFAULT_MAP),
        metadata={"description": "File with external_db mapping"}
    )


def make_mysql_url(host: str, user: str, database: str, port: str = 0, password: str = "") -> str:
    user_pass = user
    host_port = host
    if password:
        user_pass = f"{user}:{password}"
    if port:
        host_port = f"{host}:{port}"
    db_url = f"mysql://{user_pass}@{host_port}/{database}"
    return db_url


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    args = mod.args

    host = mod.args["host"]
    port = mod.args["port"]
    user = mod.args["user"]
    password = mod.args.get("password")
    database = mod.args.get("database")
    db_url = make_mysql_url(
        host=host,
        port=port,
        user=user,
        password=password,
        database=database,
    )
    dbc = DBConnection(db_url)

    external_map_path = Path(mod.args.get("external_db_map"))
    external_map = get_external_db_map(external_map_path)

    with dbc.session_scope() as session:
        seq_regions = get_seq_regions(session, external_map)

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(seq_regions, indent=2))
    else:
        print(seq_regions)


if __name__ == "__main__":
    main()
