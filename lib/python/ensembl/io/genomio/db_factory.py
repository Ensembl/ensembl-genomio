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

    if mod.args.get("output_json"):
        output_file = Path(mod.args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(dbs, indent=2))
    else:
        print(dbs)


if __name__ == "__main__":
    main()
