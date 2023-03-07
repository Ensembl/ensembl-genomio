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

from ensembl.brc4.runnable.core_server import CoreServer
from ensembl.io.genomio.features.seq_region import SeqRegion, SeqRegionAttribute, SeqRegionSynonym


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

    brc_mode = argschema.fields.Boolean(required=False, metadata={
        "description": "BRC4 mode: use organism_abbrev for species, component for division"
    })


def get_all_seq_regions(server: CoreServer, database: str) -> List[SeqRegion]:

    seq_regions = get_seq_regions(server, database)
    seqr_attribs: dict = get_seq_region_attribs(server, database)
    seqr_syns: dict = get_seq_region_synonyms(server, database)

    print(f"Got {len(seq_regions)} seq_regions")
    print(f"Got {len(seqr_attribs)} seqr_attribs")
    print(f"Got {len(seqr_syns)} seqr_syns")

    for seqr in seq_regions:
        seq_id = seqr.get_seq_region_id()
        attribs = seqr_attribs.get(seq_id)
        if attribs:
            seqr.attributes += attribs
        syns = seqr_syns.get(seq_id)
        if syns:
            seqr.synonyms += syns

    return seq_regions


def get_seq_regions(server: CoreServer, database: str) -> List[SeqRegion]:
    server.set_database(database)

    seqr_data = server.get_table_data(
        table='seq_region',
        fields=['seq_region_id', 'name', 'length', 'coord_system_id'],
    )

    seq_regions = []
    for seqr_row in seqr_data:
        seqr = SeqRegion(seq_region_id=seqr_row["seq_region_id"], name=seqr_row["name"], length=seqr_row["length"])
        seq_regions.append(seqr)
    return seq_regions


def get_seq_region_attribs(server: CoreServer, database: str) -> Dict[str, List[SeqRegionAttribute]]:
    server.set_database(database)

    seqra_data = server.get_table_data(
        table='seq_region_attrib LEFT JOIN attrib_type USING(attrib_type_id)',
        fields=['seq_region_id', 'code', 'value'],
    )

    seqr_attribs = dict()
    for row in seqra_data:
        seqr_attrib = SeqRegionAttribute(value=row["value"], code=row["code"])
        seq_id = row["seq_region_id"]
        if seq_id in seqr_attribs:
            seqr_attribs[seq_id].append(seqr_attrib)
        else:
            seqr_attribs[seq_id] = [seqr_attrib]
    return seqr_attribs


def get_seq_region_synonyms(server: CoreServer, database: str) -> Dict[str, List[SeqRegionSynonym]]:
    server.set_database(database)

    seqra_data = server.get_table_data(
        table='seq_region_synonym LEFT JOIN external_db USING(external_db_id)',
        fields=['seq_region_id', 'synonym', 'db_name'],
        constraints="synonym IS NOT NULL"
    )

    seqr_syns = dict()
    for row in seqra_data:
        seqr_syn = SeqRegionSynonym(synonym=row["synonym"], source=row["db_name"])
        seq_id = row["seq_region_id"]
        if seq_id in seqr_syns:
            seqr_syns[seq_id].append(seqr_syn)
        else:
            seqr_syns[seq_id] = [seqr_syn]
    return seqr_syns


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    args = mod.args

    server = CoreServer(
        host=mod.args["host"],
        port=args["port"],
        user=args["user"],
        password=args.get("password")
    )
    seq_regions = get_all_seq_regions(server, args["database"])
    if mod.args.get("brc_mode", 0):
        seq_regions_struct = [seqr.to_brc_dict() for seqr in seq_regions]
    else:
        seq_regions_struct = [asdict(seqr) for seqr in seq_regions]

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(seq_regions_struct, indent=2))
    else:
        print(seq_regions_struct)


if __name__ == "__main__":
    main()
