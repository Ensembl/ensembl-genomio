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
"""Fetch all Uniprot IDs mapped to proteins in a core database."""

__all__ = []

import json
import logging
from pathlib import Path
from typing import Any, Dict, List

from sqlalchemy import select
from sqlalchemy.orm import Session, joinedload

from ensembl.core.models import Translation, ObjectXref, Xref, ExternalDb
from ensembl.database import DBConnection
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def get_uniprot_map(session: Session) -> List[str]:
    """Returns a mapping from all translations to Uniprot IDs.

    Args:
        session (Session): Session from the current core.

    """
    uniprot_dbs = ("Uniprot/SPTREMBL", "Uniprot/SWISSPROT")
    xref_stmt = (
        select(Translation.stable_id, Xref.display_label, ExternalDb.db_name)
        .join(
            ObjectXref,
            (Translation.translation_id == ObjectXref.object_xref_id)
            and (ObjectXref.ensembl_object_type == "Translation"),
        )
        .join(Xref)
        .join(ExternalDb)
        .where((ExternalDb.external_db_id == Xref.external_db_id) and ExternalDb.db_name.in_(uniprot_dbs))
    )

    uniprot_maps = []
    for row in session.execute(xref_stmt).unique().all():
        tran: Translation = row[0]
        logging.warning(f"row = {row}")
        logging.warning(f"tran = {tran}")
        raise Exception()

    return uniprot_maps


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description="Fetch all the Uniprot IDs linked to translation IDs.")
    parser.add_server_arguments(include_database=True)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    dbc = DBConnection(args.url)

    with dbc.session_scope() as session:
        uniprot_map = get_uniprot_map(session)

    print(json.dumps(uniprot_map, indent=2, sort_keys=True))
