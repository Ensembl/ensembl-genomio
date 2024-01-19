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
"""Rename seq_region BRC4 names in a given core database."""

from dataclasses import dataclass
import logging
from enum import Enum, auto
from pathlib import Path
from typing import Any, Dict, List

from sqlalchemy import select, or_
from sqlalchemy.orm import Session, joinedload

from ensembl.core.models import SeqRegion, SeqRegionSynonym, SeqRegionAttrib
from ensembl.database import DBConnection
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

class Operation(Enum):
    do_nothing = auto()
    update = auto()
    insert = auto()
    not_found = auto()

    def __str__(self) -> str:
        return f"{self.name}"

@dataclass
class SeqRegionReplacement:
    """Simple record of what to rename and how."""
    name: str
    brc_name: str
    operation: Operation

    def __repr__(self) -> str:
        return f"{self.name} -> {self.brc_name} ({self.operation})"


def get_rename_map(map_file: Path) -> List[SeqRegionReplacement]:
    """Load requested renaming from a tab file."""

    seq_regions: List[SeqRegionReplacement] = []
    with map_file.open("r") as input_fh:
        for line in input_fh:
            try:
                (name, brc_name) = [item.strip() for item in line.split("\t")]
            except ValueError:
                continue
            seq_regions.append(SeqRegionReplacement(name=name, brc_name=brc_name, operation=Operation.do_nothing))
    return seq_regions


def get_seq_regions_to_replace(session: Session, rename_map: List[SeqRegionReplacement]) -> List[SeqRegionReplacement]:
    """Check each seq_region replacement against the database.

    Args:
        session: Session from the current core.
        rename_map: List of remappings to check.
    
    Returns: same list with the operation changed to the operation it can do.

    """

    for seqr in rename_map:
        seqr_stmt = (
            select(SeqRegion)
            .where(or_(SeqRegion.name == seqr.name, SeqRegionSynonym.synonym == seqr.name))
            .options(
                joinedload(SeqRegion.seq_region_synonym).joinedload(SeqRegionSynonym.external_db),
                joinedload(SeqRegion.seq_region_attrib).joinedload(SeqRegionAttrib.attrib_type),
            )
        )
        for row in session.execute(seqr_stmt).unique().all():
            db_seqr: SeqRegion = row[0]
            if not db_seqr:
                seqr.operation = Operation.not_found
            attribs = get_attribs_dict(db_seqr)
            db_brc_name = attribs.get("BRC4_seq_region_name")
            if db_brc_name:
                if seqr.brc_name == db_brc_name:
                    logging.info(f"Seq region {seqr.name} already exists with same name")
                    seqr.operation = Operation.do_nothing
                else:
                    logging.info(f"Seq region {seqr.name} already exists with a different name ({db_brc_name} instead of {seqr.brc_name})")
                    seqr.operation = Operation.update
            else:
                logging.info(f"Seq region {seqr.name} doesn't have a BRC name")
                seqr.operation = Operation.insert

    return rename_map


def get_attribs(seq_region: SeqRegion) -> List:
    """Given a seq_region, extract the attribs as value-source items.

    Args:
        seq_region (SeqRegion): The seq_region from which the attribs are extracted.

    Returns:
        All attributes as a list of dictionaries with 'value' and 'source' keys.
    """
    attribs = seq_region.seq_region_attrib
    atts = []
    if attribs:
        for attrib in attribs:
            att_obj = {"value": attrib.value, "source": attrib.attrib_type.code}
            atts.append(att_obj)
    return atts


def get_attribs_dict(seq_region: SeqRegion) -> Dict[str, Any]:
    """Extracts all the attributes of the given sequence region.

    Args:
        seq_region: Sequence region.

    Returns:
        Pairs of source and value for each attribute.
    """

    attribs = get_attribs(seq_region)
    attrib_dict = {attrib["source"]: attrib["value"] for attrib in attribs}
    return attrib_dict


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(
        description="Replace seq_region names in a core database."
    )
    parser.add_argument_src_path("input_map", help="List of seq_region names and their BRC4 name in tab format.")
    parser.add_server_arguments(include_database=True)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    rename_map = get_rename_map(args.input_map)
    print(rename_map)
    dbc = DBConnection(args.url)

    with dbc.session_scope() as session:
        seq_regions = get_seq_regions_to_replace(session, rename_map)
        print(seq_regions)
        # if args.update:
        #     logging.info("Replacing all seq_region names.")
        #     for seqr in seq_regions:
        #         logging.debug(f"Rename {seqr}")
        # else:
        #     logging.info("No change made (use --update to update the database).")


if __name__ == "__main__":
    main()