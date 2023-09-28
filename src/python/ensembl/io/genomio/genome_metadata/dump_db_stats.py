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
from sqlalchemy.orm import Session

from ensembl.database import DBConnection
from ensembl.core.models import SeqRegionAttrib, AttribType, Gene


class StatsGenerator:
    """Interface to extract stats from a core database."""

    def __init__(self, session: Session) -> None:
        self.session = session

    def get_assembly_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the assembly."""
        stats = {
            "coord_system": self.get_coord_system_tags(),
        }
        return stats

    def get_coord_system_tags(self) -> Dict[str, int]:
        """Returns a dict of stats about the coordinate systems (number of chromosomes, etc.)."""
        session = self.session

        # Assuming stats are in seq_region_attribs
        seqs_st = (
            select(SeqRegionAttrib.value, func.count(SeqRegionAttrib.value))
            .join(AttribType)
            .filter(AttribType.code == "coord_system_tag")
            .group_by(SeqRegionAttrib.value)
        )

        coords = {}
        for row in session.execute(seqs_st):
            (coord_tag, count) = row
            coords[coord_tag] = count

        return coords

    def get_annotation_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the coordinate systems (number of biotypes, etc.)."""

        stats = {
            "biotypes": self.get_biotypes(),
            "genes": self.get_genes_stats(),
        }

        return stats

    def get_biotypes(self) -> Dict[str, int]:
        """Returns a dict of stats about the genes biotypes."""
        session = self.session

        seqs_st = (
            select(Gene.biotype, func.count())
            .group_by(Gene.biotype)
        )

        biotypes = {}
        for row in session.execute(seqs_st):
            (biotype, count) = row
            biotypes[biotype] = count

        return biotypes

    def get_genes_stats(self) -> Dict[str, int]:
        """Returns a dict of stats about genes."""
        session = self.session

        totals_st = select(func.count(Gene.gene_id))
        (total,) = session.execute(totals_st).one()
        no_desc_st = select(func.count(Gene.gene_id)).filter(Gene.description == None)
        (no_desc,) = session.execute(no_desc_st).one()
        xref_desc_st = select(func.count(Gene.gene_id)).where(Gene.description.like("%[Source:%"))
        (xref_desc,) = session.execute(xref_desc_st).one()

        left_over = total - no_desc - xref_desc

        gene_stats = {
            "total": total,
            "description": {
                "empty": no_desc,
                "source_xref": xref_desc,
                "normal": left_over,
            }
        }
        return gene_stats

    def get_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the assembly and annotation."""
        all_stats = {
            "assembly": self.get_assembly_stats(),
            "annotation": self.get_annotation_stats(),
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
        generator = StatsGenerator(session)
        all_stats = generator.get_stats()

    if args.get("output_json"):
        output_file = Path(args.get("output_json"))
        with output_file.open("w") as output_fh:
            output_fh.write(json.dumps(all_stats, indent=2, sort_keys=True))
    else:
        print(json.dumps(all_stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
