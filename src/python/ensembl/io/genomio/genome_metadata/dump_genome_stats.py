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
"""Generates a JSON file representing various stats for the assembly and annotation from a core db."""

import json
from pathlib import Path
from typing import Any, Dict

from sqlalchemy import select, func
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from ensembl.core.models import SeqRegionAttrib, AttribType, Gene, Transcript
from ensembl.database import DBConnection
from ensembl.utils.argparse import ArgumentParser


class StatsGenerator:
    """Interface to extract stats from a core database."""

    def __init__(self, session: Session) -> None:
        self.session = session

    def get_assembly_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the assembly."""
        stats = {
            "coord_system": self.get_attrib_counts("coord_system_tag"),
            "locations": self.get_attrib_counts("sequence_location"),
            "codon_table": self.get_attrib_counts("codon_table"),
        }

        # Special: rename supercontigs to scaffolds for homogeneity
        stats = self._fix_scaffolds(stats)
        return stats

    def _fix_scaffolds(self, stats: Dict[str, Any]) -> Dict[str, Any]:
        """Rename supercontigs to scaffolds in a stats dict and return it."""
        coords = stats.get("coord_system", {})
        if "supercontig" in coords:
            if "scaffold" not in coords:
                coords["scaffold"] = coords["supercontig"]
                del coords["supercontig"]
        return stats

    def get_attrib_counts(self, code: str) -> Dict[str, Any]:
        """Returns a dict of count for each value counted with the attrib_type code provided.

        Args:
            code: Ensembl database attrib_type code.
        """
        session = self.session

        seqs_st = (
            select(SeqRegionAttrib.value, func.count(SeqRegionAttrib.value))
            .join(AttribType)
            .filter(AttribType.code == code)
            .group_by(SeqRegionAttrib.value)
        )

        attribs = {}
        for row in session.execute(seqs_st):
            (attrib_name, count) = row
            attribs[attrib_name] = count

        return attribs

    def get_annotation_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the coordinate systems (number of biotypes, etc.)."""

        stats = {
            "genes": self.get_feature_stats(Gene),
            "transcripts": self.get_feature_stats(Transcript),
        }

        return stats

    def get_biotypes(self, table) -> Dict[str, int]:
        """Returns a dict of stats about the feature biotypes."""
        session = self.session

        seqs_st = select(table.biotype, func.count()).group_by(table.biotype)

        biotypes = {}
        for row in session.execute(seqs_st):
            (biotype, count) = row
            biotypes[biotype] = count

        return biotypes

    def get_feature_stats(self, table) -> Dict[str, int]:
        """Returns a dict of stats about a given feature."""
        session = self.session

        totals_st = select(func.count(table.stable_id))
        (total,) = session.execute(totals_st).one()
        no_desc_st = select(func.count(table.stable_id)).filter(table.description is None)
        (no_desc,) = session.execute(no_desc_st).one()
        xref_desc_st = select(func.count(table.stable_id)).where(table.description.like("%[Source:%"))
        (xref_desc,) = session.execute(xref_desc_st).one()

        left_over = total - no_desc - xref_desc

        feat_stats = {
            "total": total,
            "biotypes": self.get_biotypes(table),
            "description": {
                "empty": no_desc,
                "source_xref": xref_desc,
                "normal": left_over,
            },
        }
        return feat_stats

    def get_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the assembly and annotation."""
        all_stats = {
            "assembly_stats": self.get_assembly_stats(),
            "annotation_stats": self.get_annotation_stats(),
        }
        return all_stats


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Fetch all the sequence regions from a core database and print them in JSON format."
    )
    parser.add_database_arguments()
    args = parser.parse_args()

    db_url = URL.create("mysql", args.user, args.password, args.host, args.port, args.database)
    dbc = DBConnection(db_url)

    with dbc.session_scope() as session:
        generator = StatsGenerator(session)
        all_stats = generator.get_stats()

    print(json.dumps(all_stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
