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
import json
import logging
from enum import Enum, auto
from pathlib import Path
from typing import Any, Dict, List

from sqlalchemy import select
from sqlalchemy.orm import Session, joinedload

from ensembl.core.models import CoordSystem, SeqRegion, SeqRegionSynonym, SeqRegionAttrib
from ensembl.database import DBConnection
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

class Operation(Enum):
    do_nothing = auto()
    update = auto()
    insert = auto()

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


def get_seq_regions(session: Session) -> List[SeqRegion]:
    """Returns all the sequence regions from the current core database.

    Include synonyms, attribs and karyotypes. Only the top level sequences are exported.

    Args:
        session (Session): Session from the current core.
        external_db_map (dict): Mapping of external_db names for the synonyms.

    """
    seq_regions = []

    seqr_stmt = (
        select(SeqRegion)
        .options(
            joinedload(SeqRegion.seq_region_synonym).joinedload(SeqRegionSynonym.external_db),
            joinedload(SeqRegion.seq_region_attrib).joinedload(SeqRegionAttrib.attrib_type),
        )
    )
    for row in session.execute(seqr_stmt).unique().all():
        seqr: SeqRegion = row[0]
        seq_region: Dict[str, Any] = {}
        seq_region = {"name": seqr.name}
        synonyms = get_synonyms(seqr)
        if synonyms:
            seq_region["synonyms"] = synonyms

        attribs = get_attribs_dict(seqr)
        if attribs:
            if "toplevel" not in attribs:
                continue
            add_attribs(seq_region, attribs)
        else:
            # Skip seq_region without attribs, not toplevel
            continue

        seq_regions.append(seq_region)

    return seq_regions


def add_attribs(seq_region: Dict, attrib_dict: Dict) -> None:
    """Map seq_regions attribs to a specific name and type defined below.

    Args:
        seq_region (Dict): A seq_region dict to modify.
        attrib_dict (Dict): The attribs for this seq_region.
    """
    bool_attribs = {
        "circular_seq": "circular",
        "non_ref": "non_ref",
    }
    int_attribs = {
        "codon_table": "codon_table",
    }
    string_attribs = {
        "BRC4_seq_region_name": "BRC4_seq_region_name",
        "EBI_seq_region_name": "EBI_seq_region_name",
        "coord_system_tag": "coord_system_level",
        "sequence_location": "location",
    }

    for name, key in bool_attribs.items():
        value = attrib_dict.get(name)
        if value:
            seq_region[key] = bool(value)

    for name, key in int_attribs.items():
        value = attrib_dict.get(name)
        if value:
            seq_region[key] = int(value)

    for name, key in string_attribs.items():
        value = attrib_dict.get(name)
        if value:
            seq_region[key] = str(value)


def get_synonyms(seq_region: SeqRegion) -> List:
    """Get all synonyms for a given seq_region. Use the mapping for synonym source names.

    Args:
        seq_region (SeqRegion): Seq_region from which the synonyms are extracted.

    Returns:
        List: All synonyms as a dict with 'name' and 'source' keys.
    """
    synonyms = seq_region.seq_region_synonym
    syns = []
    if synonyms:
        for syn in synonyms:
            if syn.external_db:
                source = syn.external_db.db_name
                syn_obj = {"name": syn.synonym, "source": source}
            else:
                syn_obj = {"name": syn.synonym}
            syns.append(syn_obj)

    return syns


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

def get_rename_map(map_file: Path) -> List[SeqRegionReplacement]:
    """Load requested renaming from a tab file."""

    seq_regions: List[SeqRegionReplacement] = []
    with map_file.open("r") as input_fh:
        for line in input_fh:
            (name, brc_name) = [item.strip() for item in line.split("\t")]
            seq_regions.append(SeqRegionReplacement(name=name, brc_name=brc_name, operation=Operation.do_nothing))
    return seq_regions

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
    # dbc = DBConnection(args.url)

    # with dbc.session_scope() as session:
        # seq_regions = get_seq_regions_to_replace(session, rename_map)
        # if args.update:
        #     logging.info("Replacing all seq_region names.")
        #     for seqr in seq_regions:
        #         logging.debug(f"Rename {seqr}")
        # else:
        #     logging.info("No change made (use --update to update the database).")


if __name__ == "__main__":
    main()