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
"""Update descriptions from a functional annotation file into a core database."""

__all__ = [
    "get_core_data",
    "load_descriptions",
]

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from sqlalchemy.orm import Session
from sqlalchemy import and_, select

import ensembl.io.genomio
from ensembl.core.models import Gene, Transcript, ObjectXref, Xref
from ensembl.io.genomio.utils import get_json
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.database import DBConnection
from ensembl.utils.logging import init_logging_with_args


FEAT_TABLE = {
    "gene": "gene",
    "mobile_element": "gene",
    "transcript": "transcript",
}

FeatStruct = Tuple[str, str, str]


def get_core_data(session: Session, table: str, match_xrefs: bool = False) -> Dict[str, FeatStruct]:
    """Returns the table descriptions from a core database.

    Args:
        session: Session open on a core database.
        table: "gene" or "transcript" table from the core database.
        match_xrefs: If the IDs do not match, try to match an Xref ID instead.
    """

    if table == "gene":
        stmt = (
            select(Gene.gene_id, Gene.stable_id, Gene.description, Xref.dbprimary_acc)
            .select_from(Gene)
            .outerjoin(
                ObjectXref,
                and_(Gene.gene_id == ObjectXref.ensembl_id, ObjectXref.ensembl_object_type == "gene"),
            )
            .outerjoin(Xref)
        )
    elif table == "transcript":
        stmt = (
            select(Transcript.transcript_id, Transcript.stable_id, Transcript.description, Xref.dbprimary_acc)
            .select_from(Transcript)
            .outerjoin(
                ObjectXref,
                and_(
                    Transcript.transcript_id == ObjectXref.ensembl_id,
                    ObjectXref.ensembl_object_type == "transcript",
                ),
            )
            .outerjoin(Xref)
        )
    else:
        raise ValueError(f"Table {table} is not supported")

    feat_data = {}
    for row in session.execute(stmt):
        (feat_id, stable_id, desc, xref_name) = row
        feat_struct: FeatStruct = (feat_id, stable_id, desc)
        feat_data[stable_id.lower()] = feat_struct
        if match_xrefs and xref_name:
            feat_data[xref_name.lower()] = feat_struct

    return feat_data


def load_descriptions(
    session: Session,
    func_file: Path,
    report: bool = False,
    do_update: bool = False,
    match_xrefs: bool = True,
) -> None:
    """Loads gene and transcript descriptions into a core database.

    Args:
        session: Session open on a core database.
        func_file: JSON file with the annotation information.
        report: Print the mapping of changes to perform in the standard output.
        do_update: Actually update the core database.
        match_xrefs: If the IDs do not match, try to match an Xref ID instead.
    """
    func = get_json(func_file)
    logging.info(f"{len(func)} annotations from {func_file}")
    table_to_update = {"gene": Gene, "transcript": Transcript}
    for table, mapped_table in table_to_update.items():
        logging.info(f"Checking {table} descriptions")
        feat_func = [feat for feat in func if feat["object_type"] == table]
        logging.info(f"{len(feat_func)} {table} annotations from {func_file}")
        feat_data = get_core_data(session, table, match_xrefs)
        logging.info(f"Loaded {len(feat_data)} {table} data")

        stats = {
            "not_supported": 0,
            "not_found": 0,
            "same": 0,
            "same_empty": 0,
            "empty_but_xref": 0,
            "to_update_replace": 0,
            "to_update_remove": 0,
        }
        # Compare, only keep the descriptions that have changed
        features_to_update = _get_features_to_update(
            table, feat_func, feat_data, stats, report=report, do_update=do_update, match_xrefs=match_xrefs
        )

        # Show stats for this feature type
        for stat, count in stats.items():
            if count == 0:
                continue
            logging.info(f"{stat} = {count}")

        if do_update:
            logging.info(f"Now updating {len(features_to_update)} rows...")
            session.bulk_update_mappings(mapped_table, features_to_update)
            session.commit()


def _get_cur_feat(
    feat_data: Dict[str, FeatStruct], new_feat: Dict[str, Any], match_xrefs: bool = False
) -> Optional[FeatStruct]:
    """Match a feature ID, synonyms or xrefs to a core stable ID and return the matching core feature.

    Returns None if no match.
    """
    # Match with the ID
    cur_feat = feat_data.get(new_feat["id"].lower())

    # Fall back to a synonym
    if not cur_feat and "synonyms" in new_feat:
        for syn in new_feat["synonyms"]:
            cur_feat = feat_data.get(syn.lower())
            if cur_feat:
                break

    # Fall back to an xref
    if not cur_feat and match_xrefs and "xrefs" in new_feat:
        for xref in new_feat["xrefs"]:
            cur_feat = feat_data.get(xref["id"].lower())
            if cur_feat:
                break

    return cur_feat


def _get_features_to_update(
    table: str,
    feat_func: List[Dict[str, Any]],
    feat_data: Dict[str, FeatStruct],
    stats: Dict[str, int],
    *,
    report: bool = False,
    do_update: bool = False,
    match_xrefs: bool = True,
) -> List[Dict[str, Any]]:
    """Checks a list of features and returns those whose description we want to update.

    Args:
        table: "gene" or "transcript" table for the features.
        feat_func: The features to check.
        feat_data: The features in the database.
        stats: Record the number of features checked in different cases.
        report: Print a report line for each feature to standard output.
        do_update: Actually update the database.
        match_xrefs: Use xref IDs if feature ID does not match a feature in the database.

    Returns:
        The list of features with their operation changed to update or insert.
    """
    to_update = []
    for new_feat in feat_func:
        cur_feat = _get_cur_feat(feat_data, new_feat, match_xrefs)

        # No match in the end
        if not cur_feat:
            logging.debug(f"Not found: {table} '{new_feat['id']}'")
            stats["not_found"] += 1
            continue

        # Prepare some data to compare
        new_stable_id = new_feat["id"]
        new_desc = new_feat.get("description", "")
        (row_id, cur_stable_id, cur_desc) = cur_feat

        # No description: replace unless the current description is from an Xref
        if not cur_desc:
            cur_desc = ""
        if not new_desc:
            if cur_desc == "":
                stats["same_empty"] += 1
                continue
            if "[Source:" in cur_desc:
                stats["empty_but_xref"] += 1
                continue
            stats["to_update_remove"] += 1

        # Compare the descriptions
        elif new_desc == cur_desc:
            stats["same"] += 1
            continue
        # At this point, we have a new description to update
        else:
            stats["to_update_replace"] += 1

        # Directly print the mapping
        if report:
            line = (table, new_stable_id, cur_stable_id, cur_desc, new_desc)
            print("\t".join(line))

        # Add to the batch list of updates for the core db
        if do_update:
            update_key = f"{table}_id"
            to_update.append({update_key: row_id, "description": new_desc})

    return to_update


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_server_arguments(include_database=True)
    parser.add_argument_src_path("--func_file", required=True, help="Input functional annotation JSON")
    parser.add_argument("--report", action="store_true", help="Show what change would be made")
    parser.add_argument("--update", action="store_true", help="Make the changes to the database")
    parser.add_argument(
        "--match_xrefs", action="store_true", help="Use xref IDs to match features if IDs do not work"
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnection(args.url)
    with dbc.session_scope() as session:
        load_descriptions(session, args.func_file, args.report, args.update, match_xrefs=args.match_xrefs)
