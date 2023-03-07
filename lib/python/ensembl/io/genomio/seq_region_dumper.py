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
from ensembl.io.genomio.features.seq_region import CoordSystem, SeqRegion, SeqRegionAttribute, SeqRegionSynonym, KaryotypeBand


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

    coord_systems = get_coord_systems(server, database)
    print(f"Got {len(coord_systems)} coord_systems")
    seq_regions = get_seq_regions(server, database, coord_systems)
    print(f"Got {len(seq_regions)} seq_regions")
    seqr_attribs: dict = get_seq_region_attribs(server, database)
    print(f"Got {len(seqr_attribs)} seqr_attribs")
    seqr_syns: dict = get_seq_region_synonyms(server, database)
    print(f"Got {len(seqr_syns)} seqr_syns")
    karyotypes: dict = get_karyotypes(server, database)
    print(f"Got {len(karyotypes)} karyotypes")

    final_seq_regions = []
    for seqr in seq_regions:
        # Add Coord system id
        coord_id = seqr.coord_system_id
        coord = coord_systems.get(coord_id)
        if coord:
            seqr.coord_system = coord

        # Add attribs
        seq_id = seqr.get_seq_region_id()
        attribs = seqr_attribs.get(seq_id)
        if attribs:
            seqr.attributes += attribs

        # Add Synonyms
        syns = seqr_syns.get(seq_id)
        if syns:
            seqr.synonyms += syns

        # Add Karyotypes
        kars = karyotypes.get(seq_id)
        if kars:
            seqr.karyotype_bands += kars

        # Filtering
        if seqr.is_top_level():
            final_seq_regions.append(seqr)

    print(f"Got {len(final_seq_regions)} seq_regions dumped")

    return final_seq_regions


def get_coord_systems(server: CoreServer, database: str) -> Dict[str, CoordSystem]:
    server.set_database(database)

    seqr_data = server.get_table_data(
        table='coord_system',
        fields=['coord_system_id', 'species_id', 'name', 'version', 'attrib'],
    )

    coord_systems = {}
    for seqr_row in seqr_data:
        attrib = seqr_row.get("attrib")
        if not attrib:
            attrib = []
        else:
            attrib = list(attrib)
        if not seqr_row["version"]:
            seqr_row["version"] = ""
        if "default_version" in attrib:
            coord = CoordSystem(
                coord_system_id=seqr_row["coord_system_id"],
                species_id=seqr_row["species_id"],
                name=seqr_row["name"],
                version=seqr_row["version"],
                attrib=attrib,
            )
            coord_systems[coord.coord_system_id] = coord
    return coord_systems


def get_seq_regions(
        server: CoreServer, database: str, coord_systems: Dict[str, CoordSystem]) -> List[SeqRegion]:
    server.set_database(database)

    coords = ", ".join([str(c.coord_system_id) for c in coord_systems.values()])
    if coords:
        constraints = f"coord_system_id IN ({coords})"

    seqr_data = server.get_table_data(
        table='seq_region',
        fields=['seq_region_id', 'name', 'length', 'coord_system_id'],
        constraints=constraints
    )

    seq_regions = []
    for seqr_row in seqr_data:
        seqr = SeqRegion(
            seq_region_id=seqr_row["seq_region_id"],
            name=seqr_row["name"],
            length=seqr_row["length"],
            coord_system_id=seqr_row["coord_system_id"]
        )
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


def get_karyotypes(server: CoreServer, database: str) -> Dict[str, List[KaryotypeBand]]:
    server.set_database(database)

    kar_data = server.get_table_data(
        table='karyotype',
        fields=['seq_region_id', 'seq_region_start', 'seq_region_end', 'band', 'stain'],
    )

    karyotype = dict()
    for row in kar_data:
        kar = KaryotypeBand(
            band=row["band"],
            seq_region_start=row["seq_region_start"],
            seq_region_end=row["seq_region_end"],
            stain=row["stain"],
        )
        seq_id = row["seq_region_id"]
        if seq_id in karyotype:
            karyotype[seq_id].append(kar)
        else:
            karyotype[seq_id] = [kar]
    return karyotype


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
