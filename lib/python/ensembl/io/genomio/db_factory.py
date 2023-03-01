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

from ensembl.brc4.runnable.core_server import CoreServer


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    # Server parameters
    host = argschema.fields.String(metadata={
        "required": True, "description": "Host to the server with EnsEMBL databases"
    })
    port = argschema.fields.Integer(metadata={
        "required": True, "description": "Port to use"
    })
    host = argschema.fields.String(metadata={
        "required": True, "description": "Host to use"
    })
    user = argschema.fields.String(metadata={
        "required": True, "description": "User to use"
    })
    password = argschema.fields.String(metadata={
        "required": False, "description": "Password to use"
    })

    # Filters
    prefix = argschema.fields.String(metadata={
        "required": False, "description": "Prefix to filter the databases"
    })
    build = argschema.fields.String(metadata={
        "required": False, "description": "Build to filter the databases"
    })
    version = argschema.fields.String(metadata={
        "required": False, "description": "EnsEMBL version to filter the databases"
    })
    brc_mode = argschema.fields.Boolean(metadata={
        "required": False, "description": "BRC4 mode: use organism_abbrev for species, component for division"
    })


def format_db_data(server: CoreServer, dbs: List[str], brc_mode: False) -> List[Dict]:
    db_datas = list()
    for db in dbs:
        server.set_database(db)
        metadata = server.get_db_metadata()

        species = get_metadata_value(metadata, "species.production_name")
        division = get_metadata_value(metadata, "species.division")

        if brc_mode:
            brc_organism = get_metadata_value(metadata, "BRC4.organism_abbrev")
            brc_component = get_metadata_value(metadata, "BRC4.component")
            if brc_organism:
                species = brc_organism
            if brc_component:
                division = brc_component
        
        if not division:
            division = 'all'
        
        db_data = {
            "database": db,
            "species": species,
            "division": division,
        }
        db_datas.append(db_data)
    return db_datas


def get_metadata_value(metadata, key) -> str:
    if key in metadata and len(metadata[key]) > 0:
        return metadata[key][0]
    else:
        return None


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)

    server = CoreServer(
        host=mod.args["host"],
        port=mod.args["port"],
        user=mod.args["user"],
        password=mod.args.get("password")
    )

    prefix = mod.args.get("prefix")
    build = mod.args.get("build")
    version = mod.args.get("version")
    dbs = server.get_cores(prefix=prefix, build=build, version=version)

    brc_mode = mod.args.get("brc_mode")
    dbs_data = format_db_data(server, dbs, brc_mode)

    if mod.args.get("output_json"):
        output_file = Path(mod.args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(dbs_data, indent=2))
    else:
        print(dbs)


if __name__ == "__main__":
    main()
