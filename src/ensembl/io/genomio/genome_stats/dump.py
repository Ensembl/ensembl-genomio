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
"""Generates a JSON representation of the genome stats (assembly and annotation) from a core database."""

__all__ = ["StatsGenerator"]

from dataclasses import dataclass
import json
from typing import Any, Dict

from sqlalchemy import select, func
from sqlalchemy.orm import Session

from ensembl.core.models import SeqRegionAttrib, AttribType, Gene, Transcript
import ensembl.io.genomio
from ensembl.io.genomio.database import DBConnectionLite
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.database import StrURL
from ensembl.utils.logging import init_logging_with_args


@dataclass
class StatsGenerator:
    """Interface to extract genome stats from a core database."""

    session: Session

    def get_assembly_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the assembly."""
        stats = {
            "coord_system": self.get_attrib_counts("coord_system_tag"),
            "locations": self.get_attrib_counts("sequence_location"),
            "codon_table": self.get_attrib_counts("codon_table"),
        }
        # Special: rename supercontigs to scaffolds for homogeneity
        StatsGenerator._fix_scaffolds(stats)
        return stats

    @staticmethod
    def _fix_scaffolds(stats: Dict[str, Any]) -> None:
        """Renames supercontigs to scaffolds in the provided stats.

        If scaffolds are present already, nothing is done.

        Args:
            stats: Statistics dictionary.

        """
        coords = stats.get("coord_system", {})
        if "supercontig" in coords and "scaffold" not in coords:
            coords["scaffold"] = coords["supercontig"]
            del coords["supercontig"]

    def get_attrib_counts(self, code: str) -> Dict[str, Any]:
        """Returns a dict of count for each value counted with the attrib_type code provided.

        Args:
            code: Ensembl database attrib_type code.

        """
        seqs_st = (
            select(SeqRegionAttrib.value, func.count())  # pylint: disable=not-callable
            .join(AttribType)
            .filter(AttribType.code == code)
            .group_by(SeqRegionAttrib.value)
        )
        attributes = {}
        for row in self.session.execute(seqs_st):
            (attribute_name, count) = row
            attributes[attribute_name] = count
        return attributes

    def get_annotation_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the coordinate systems (number of biotypes, etc.)."""
        stats = {
            "genes": self.get_feature_stats(Gene),
            "transcripts": self.get_feature_stats(Transcript),
        }
        return stats

    def get_biotypes(self, table: Any) -> Dict[str, int]:
        """Returns a dict of stats about the feature biotypes."""
        # pylint: disable-next=not-callable
        seqs_st = select(table.biotype, func.count()).group_by(table.biotype)
        biotypes = {}
        for row in self.session.execute(seqs_st):
            (biotype, count) = row
            biotypes[biotype] = count
        return biotypes

    def get_feature_stats(self, table: Any) -> Dict[str, int]:
        """Returns a dict of stats about a given feature."""
        session = self.session
        totals_st = select(func.count()).select_from(table)  # pylint: disable=not-callable
        (total,) = session.execute(totals_st).one()
        # pylint: disable-next=singleton-comparison,not-callable
        no_desc_st = select(func.count()).filter(table.description.is_(None))
        (no_desc,) = session.execute(no_desc_st).one()
        # pylint: disable-next=not-callable
        xref_desc_st = select(func.count()).where(table.description.like("%[Source:%"))
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

    def get_genome_stats(self) -> Dict[str, Any]:
        """Returns a dict of stats about the assembly and annotation."""
        genome_stats = {
            "assembly_stats": self.get_assembly_stats(),
            "annotation_stats": self.get_annotation_stats(),
        }
        return genome_stats


def dump_genome_stats(url: StrURL) -> Dict[str, Any]:
    """Returns JSON object containing the genome stats (assembly and annotation) of the given core database.

    Args:
        url: Core database URL.

    """
    dbc = DBConnectionLite(url)
    with dbc.session_scope() as session:
        generator = StatsGenerator(session)
        genome_stats = generator.get_genome_stats()
        return genome_stats


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_server_arguments(include_database=True)
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    genome_stats = dump_genome_stats(args.url)
    print(json.dumps(genome_stats, indent=2, sort_keys=True))
