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
"""Load functional annotation from a file into a core database.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List

from ensembl.database import DBConnection
from ensembl.core.models import Gene, Transcript, ObjectXref, Xref
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args
from sqlalchemy import select, and_, or_
from sqlalchemy.orm import Session

from ensembl.io.genomio.utils import get_json


FEAT_TABLE = {
    "gene": "gene",
    "mobile_element": "gene",
    "transcript": "transcript",
}


def load_func(session: Session, func_file: Path, report: bool = False, update: bool = False) -> None:
    """Load functional annotation in a core, can report or replace."""

    func = get_json(func_file)
    logging.info(f"{len(func)} annotations from {func_file}")

    stats = {
        "not_supported": 0,
        "not_found": 0,
        "differs": 0,
        "same": 0,
        "no_new": 0,
    }
    all_num = 0
    for feat in func:
        all_num += 1
        if all_num % 500 == 0:
            logging.info(f"Processed {all_num} features")

        if "description" not in feat:
            stats["no_new"] += 1
            continue
        logging.debug(f"Look for {feat}")
        # Check which table to compare to
        try:
            table = FEAT_TABLE[feat["object_type"]]
        except KeyError:
            logging.debug(f"Ignore feat {feat} not supported for loading")
            stats["not_supported"] += 1
            continue
        
        # Get the feature from the database
        if table == "gene":
            stmt = (
                select(Gene)
                .select_from(Gene)
                .join(ObjectXref, Gene.gene_id == ObjectXref.ensembl_id)
                .join(Xref)
                .where(ObjectXref.ensembl_object_type == "gene")
                .where(or_(Gene.stable_id == feat["id"], Xref.dbprimary_acc == feat["id"]))
            )
            gene = session.scalars(stmt).one()
            if not gene:
                logging.debug(f"No gene found with this name: {feat['id']}")
                stats["not_found"] += 1
                continue

            gene_description = "(NULL)"
            if gene.description:
                gene_description = gene.description
            
            # Check that the description differs
            new_desc = feat["description"]
            if new_desc == gene_description:
                logging.debug(f"Gene {gene.stable_id} has the same description")
                stats["same"] += 1
                continue
            logging.debug(f"Gene {gene.stable_id} has a different description")
            stats["differs"] += 1
            if report:
                line = (feat['id'], gene.stable_id, gene_description, feat['description'])
                print("\t".join(line))
            if update:
                logging.debug("TODO: add actual change")

        elif table == "transcript":
            stmt = (
                select(Transcript)
                .select_from(Transcript)
                .join(ObjectXref, Transcript.transcript_id == ObjectXref.ensembl_id)
                .join(Xref)
                .where(ObjectXref.ensembl_object_type == "transcript")
                .where(or_(Transcript.stable_id == feat["id"], Xref.dbprimary_acc == feat["id"]))
            )
            transc = session.scalars(stmt).one()
            if not transc:
                logging.debug(f"No transcript found with this name: {feat['id']}")
                stats["not_found"] += 1
                continue
            
            # Check that the description differs
            new_desc = feat["description"]
            if new_desc == transc.description:
                logging.debug(f"Transcript {transc.stable_id} has the same description")
                stats["same"] += 1
                continue
            logging.debug(f"Transcript {transc.stable_id} has a different description")
            stats["differs"] += 1
            if report:
                line = (feat['id'], transc.stable_id, gene.description, feat['description'])
                print("\t".join(line))
            if update:
                logging.debug("TODO: add actual change")
        else:
            logging.info(f"Not checking {feat}: not table")
    
    for stat, count in stats.items():
        if count == 0:
            continue
        logging.info(f"{stat} = {count}")



def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_server_arguments(include_database=True)
    parser.add_argument_src_path("--func_file", required=True, help="Input Functional annotation json")
    parser.add_argument("--report", action="store_true", help="Show what change would be made")
    parser.add_argument("--update", action="store_true", help="Make the changes to the database")
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnection(args.url)
    with dbc.session_scope() as session:
        load_func(session, args.func_file, args.report, args.update)

if __name__ == "__main__":
    main()
