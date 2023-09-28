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
"""Generates a JSON file representing various stats for the assembly and annotation from a core db.
"""

import json
from pathlib import Path
from typing import Any, Dict

import argschema
from sqlalchemy import select, func
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session, joinedload

from ensembl.database import DBConnection
from ensembl.core.models import SeqRegion, SeqRegionAttrib, AttribType


def get_assembly_stats(session: Session) -> Dict[str, Any]:
    stats = {
        "coord_system": get_coord_system_tags(session),
    }
    return stats


def get_coord_system_tags(session: Session) -> Dict[str, int]:
    # Assuming stats are in seq_region_attribs
    seqs_st = (
        select(SeqRegionAttrib.value, func.count(SeqRegionAttrib.value))
        .join(AttribType)
        .filter(AttribType.code == "coord_system_tag")
        .group_by(SeqRegionAttrib.value)
    )

    coords = {}
    for row in session.execute(seqs_st):
        (coord, value) = row
        coords[coord] = value

    return coords


def get_annotation_stats(session: Session) -> Dict[str, Any]:
    return {}


def get_stats(session: Session) -> Dict[str, Any]:
    all_stats = {
        "assembly": get_assembly_stats(session),
        "annotation": get_annotation_stats(session),
    }
    return all_stats


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

    with dbc.session_scope() as session:
        all_stats = get_stats(session)

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(all_stats, indent=2, sort_keys=True))
    else:
        print(json.dumps(all_stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
