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
"""Rename seq_region BRC names in a given core database."""

__all__ = [
    "Operation",
    "SeqRegionReplacement",
    "get_rename_map",
    "get_seq_regions_to_replace",
    "rename_seq_regions",
    "update_seq_region_name",
]

from dataclasses import dataclass
import logging
from enum import Enum, auto
from pathlib import Path
from typing import Dict, List

from sqlalchemy import select, or_
from sqlalchemy.orm import Session

from ensembl.core.models import SeqRegion, SeqRegionSynonym, SeqRegionAttrib
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.database import DBConnection
from ensembl.utils.logging import init_logging_with_args


class Operation(Enum):
    """Controlled list of operations for the seq_region replacement."""

    DO_NOTHING = auto()
    UPDATE = auto()
    INSERT = auto()
    NOT_FOUND = auto()

    def __str__(self) -> str:
        return f"{self.name}"


@dataclass
class SeqRegionReplacement:
    """Simple record of what to rename and how."""

    name: str
    brc_name: str
    operation: Operation
    old_brc_name: str = ""
    seq_region_id: int = 0

    def __repr__(self) -> str:
        return f"{self.name}: '{self.old_brc_name}' -> '{self.brc_name}' ({self.operation})"


def get_rename_map(map_file: Path) -> List[SeqRegionReplacement]:
    """Load requested renaming from a tab file."""

    seq_regions: List[SeqRegionReplacement] = []
    with map_file.open("r") as input_fh:
        for line in input_fh:
            try:
                (name, brc_name) = [item.strip() for item in line.split("\t")]
            except ValueError:
                logging.warning(f"Skipping line '{line.strip()}', wrong format")
                continue
            seq_regions.append(
                SeqRegionReplacement(name=name, brc_name=brc_name, operation=Operation.DO_NOTHING)
            )
    return seq_regions


def get_seq_regions_to_replace(
    session: Session, rename_map: List[SeqRegionReplacement]
) -> List[SeqRegionReplacement]:
    """Check each seq_region replacement against the database.

    Args:
        session: Session from the current core.
        rename_map: List of remappings to check.

    Returns:
        Same list of remappings with the operation changed to the operation it can do.

    """

    seq_regions = []
    for seqr in rename_map:
        seqr_stmt = (
            select(SeqRegion)
            .select_from(SeqRegion)
            .outerjoin(SeqRegion.seq_region_attrib)
            .outerjoin(SeqRegion.seq_region_synonym)
            .where(or_(SeqRegion.name == seqr.name, SeqRegionSynonym.synonym == seqr.name))
        )

        rows = session.execute(seqr_stmt).unique().all()
        if not rows:
            logging.warning(f"Seq region not found: {seqr}")
            seqr.operation = Operation.NOT_FOUND
            continue
        for row in rows:
            db_seqr: SeqRegion = row[0]
            if not db_seqr:
                seqr.operation = Operation.NOT_FOUND
                continue
            seqr.seq_region_id = db_seqr.seq_region_id

            attribs = _get_attribs(db_seqr)
            db_brc_name = attribs.get("BRC4_seq_region_name", "")
            seqr.old_brc_name = db_brc_name
            if not db_brc_name:
                logging.info(f"Seq region '{seqr.name}' doesn't have a BRC name")
                seqr.operation = Operation.INSERT
            elif seqr.brc_name == db_brc_name:
                logging.info(f"Seq region {seqr.name} already exists with same name")
                seqr.operation = Operation.DO_NOTHING
                continue
            else:
                logging.info(
                    f"Seq region {seqr.name} already exists as {db_brc_name} instead of {seqr.brc_name}"
                )
                seqr.operation = Operation.UPDATE
            seq_regions.append(seqr)

    return seq_regions


def _get_attribs(seq_region: SeqRegion) -> Dict[str, str]:
    """Given a seq_region, extract the attribs as value-source items.

    Args:
        seq_region: The seq_region from which the attribs are extracted.

    Returns:
        Pairs of source and value for each attribute.
    """
    db_attribs = seq_region.seq_region_attrib
    attribs_dict = {}
    for attrib in db_attribs:
        attribs_dict[attrib.attrib_type.code] = attrib.value
    return attribs_dict


def update_seq_region_name(
    session: DBConnection, seq_region: SeqRegionReplacement, update: bool = False
) -> None:
    """Update BRC seq_region based on the operation provided.

    Args:
        session: SQLAlchemy session.
        seq_region: Seq region replacement.
        update: Make actual changes to the database.

    """
    if not update:
        logging.info(f"Rename {seq_region} (fake update)")
        return
    logging.info(f"Rename {seq_region}")

    brc_attrib_id = 548
    if seq_region.operation == Operation.UPDATE:
        stmt = select(SeqRegionAttrib).where(
            SeqRegionAttrib.seq_region_id == seq_region.seq_region_id,
            SeqRegionAttrib.attrib_type_id == brc_attrib_id,
        )
        seq_attrib = session.scalars(stmt).one()
        seq_attrib.value = seq_region.brc_name
    elif seq_region.operation == Operation.INSERT:
        logging.warning(f"Not supported: insertion for {seq_region}")
    else:
        logging.warning(f"Cannot update/insert seq_region without clear operation {seq_region}")


def rename_seq_regions(dbc: DBConnection, input_map: Path, do_update: bool = False) -> None:
    """Rename seq_regions in a core from a list of replacements.

    Args:
        dbc: Connection to a core database.
        input_map: Path to a file with 2 columns: current_name, new_name.
        do_update: Flag to do actual change to the database.
    """
    rename_map = get_rename_map(input_map)

    with dbc.session_scope() as session:
        seq_regions = get_seq_regions_to_replace(session, rename_map)
        logging.info(f"Replacing {len(seq_regions)} seq_region names")
        for seqr in seq_regions:
            update_seq_region_name(session, seqr, do_update)
        if not do_update:
            logging.warning("No change made (use --update to update the database).")


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path(
        "input_map", help="List of seq_region names and their BRC name in tab format."
    )
    parser.add_argument("--update", action="store_true", help="Make changes to the database.")
    parser.add_server_arguments(include_database=True)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnection(args.url)
    rename_seq_regions(dbc, args.input_map, args.update)
