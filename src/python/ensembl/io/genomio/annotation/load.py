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


def get_core_data(session: Session, table: str) -> Dict[str, str]:
    """Load descriptions from a core."""

    if table == "gene":
        stmt = (
            select(Gene.stable_id, Gene.description, Xref.dbprimary_acc)
            .select_from(Gene).join(ObjectXref, Gene.gene_id == ObjectXref.ensembl_id)
            .join(Xref)
            .where(ObjectXref.ensembl_object_type == "gene")
        )
    elif table == "transcript":
        stmt = (
            select(Transcript.stable_id, Transcript.description, Xref.dbprimary_acc)
            .select_from(Transcript).join(ObjectXref, Transcript.transcript_id == ObjectXref.ensembl_id)
            .join(Xref)
            .where(ObjectXref.ensembl_object_type == "transcript")
        )

    feat_data = {}
    for row in session.execute(stmt):
        (name, desc, xref_name) = row
        feat_struct = (name, desc)
        feat_data[name] = feat_struct
        feat_data[xref_name] = feat_struct
    return feat_data


def load_descriptions(session: Session, func_file: Path, report: bool = False, update: bool = False) -> None:
    """Load gene and transcript descriptions in a core database.
    """
    func = get_json(func_file)
    logging.info(f"{len(func)} annotations from {func_file}")
    for table in ("gene", "transcript"):
        logging.info(f"Checking {table} descriptions")
        feat_func = [feat for feat in func if feat["object_type"] == table]
        logging.info(f"{len(feat_func)} {table} annotations from {func_file}")
        feat_data = get_core_data(session, table)
        logging.info(f"Loaded {len(feat_data)} {table} data")

        stats = {
            "not_supported": 0,
            "not_found": 0,
            "differs": 0,
            "same": 0,
            "no_new": 0,
        }
        # Compare, only keep the descriptions that have changed
        for new_feat in feat_func:
            if "description" not in new_feat:
                stats["no_new"] += 1
                continue
            try:
                current_feat = feat_data[new_feat["id"]]
            except KeyError:
                stats["not_found"] += 1
                continue
            
            new_id = new_feat["id"]
            cur_id = current_feat[0]
            new_desc = new_feat["description"]
            cur_desc = current_feat[1]
            if not cur_desc:
                cur_desc = ""

            if new_desc == cur_desc:
                stats["same"] += 1
                continue
            
            stats["differs"] += 1
            if report:
                line = (table, new_id, cur_id, cur_desc, new_desc)
                print("\t".join(line))

        # Show stats for this feature type
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
        load_descriptions(session, args.func_file, args.report, args.update)

if __name__ == "__main__":
    main()
