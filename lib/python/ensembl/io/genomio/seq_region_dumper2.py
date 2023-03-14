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

from dataclasses import asdict
import json
from pathlib import Path
from typing import Dict, List

import argschema
from sqlalchemy import select
from sqlalchemy.orm import Session

from ensembl.database import DBConnection
from ensembl.core.models import CoordSystem, SeqRegion, SeqRegionSynonym, SeqRegionAttrib

ROOT_DIR = Path(__file__).parent / "../../../../.."
DEFAULT_MAP = ROOT_DIR / "data/external_db_map_default.txt"


def get_seq_regions(session: Session) -> List[SeqRegion]:
    seqr_stmt = select(SeqRegion)
    seq_regions = []
    for row in session.execute(seqr_stmt):
        seqr: SeqRegion = row[0]
        seq_region = dict()
        seq_region = {
            "name": seqr.name,
            "length": seqr.length
        }
        synonyms = seqr.synonyms
        seq_regions.append(seq_region)
    return seq_regions


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

    external_map = Path(mod.args.get("external_db_map"))

    with dbc.session_scope() as session:
        seq_regions = get_seq_regions(session)

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(seq_regions, indent=2))
    else:
        print(seq_regions)


if __name__ == "__main__":
    main()
