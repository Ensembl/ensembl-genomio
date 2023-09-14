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
"""Generates a JSON file representing the genome metadata from a core database.
"""

import json
from pathlib import Path
from typing import Any, Dict

import argschema
from sqlalchemy import select
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from ensembl.database import DBConnection
from ensembl.core.models import Meta


def get_genome_metadata(session: Session) -> Dict[str, Any]:
    """Retrieve a select list of metadata from the core database.

    Args:
        session: Session for the current core.

    Returns:
        A nested dict.
    """
    gmeta: Dict[str, Any] = {}

    gmeta_st = select(Meta)
    for row in session.execute(gmeta_st).unique().all():
        dat = row[0]
        meta_key = dat.meta_key
        meta_value = dat.meta_value

        if "." in meta_key:
            (high_key, low_key) = meta_key.split(".")
            if high_key in gmeta:
                if low_key in gmeta[high_key]:
                    gmeta[high_key][low_key].append(meta_value)
                else:
                    gmeta[high_key][low_key] = [meta_value]
            else:
                gmeta[high_key] = {}
                gmeta[high_key][low_key] = [meta_value]
        else:
            if meta_key in gmeta:
                gmeta[meta_key].append(meta_value)
            else:
                gmeta[meta_key] = [meta_value]

    return gmeta


def filter_genome_meta(gmeta: Dict[str, Any]) -> Dict[str, Any]:
    meta_list = {
        "species": {"taxonomy_id"},
        "BRC4": {"organism_abbrev", "component"},
    }

    gmeta_out: Dict[str, Any] = {}
    for mkey, mval in meta_list.items():
        if mkey in gmeta:
            gmeta_out[mkey] = gmeta[mkey]
    
    return gmeta_out




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
        genome_meta = get_genome_metadata(session)
        genome_meta = filter_genome_meta(genome_meta)

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(genome_meta, indent=2, sort_keys=True))
    else:
        print(json.dumps(genome_meta, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
