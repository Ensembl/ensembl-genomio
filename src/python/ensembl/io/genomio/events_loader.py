#!/usr/bin/env python3
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


from pathlib import Path
from typing import List
from dataclasses import dataclass


import argschema
from sqlalchemy.orm import Session

from ensembl.database import DBConnection


@dataclass
class IdEvent:
    from_id: str
    to_id: str
    release: str
    release_date: str


def load_events(input_file: Path) -> List[IdEvent]:
    events: List[IdEvent] = []

    with input_file.open('r') as events_fh:
        for line in events_fh:
            line.strip()
            if line == "":
                continue
            (from_id, to_id, _, release, release_date) = line.split("\t")
            event = IdEvent(
                from_id=from_id,
                to_id=to_id,
                release=release,
                release_date=release_date
            )
            events.append(event)
    return events


def write_events(session: Session, events: List[IdEvent]) -> None:
    return


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
    input_file = argschema.fields.InputFile(required=True, metadata={"description": "Input file"})


def make_mysql_url(host: str, user: str, database: str, port: int = 0, password: str = "") -> str:
    user_pass = user
    host_port = host
    if password:
        user_pass = f"{user}:{password}"
    if port:
        host_port = f"{host}:{port}"
    db_url = f"mysql://{user_pass}@{host_port}/{database}"
    return db_url


def main() -> None:
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)

    # Start
    host = mod.args["host"]
    port = int(mod.args["port"])
    user = mod.args["user"]
    password = mod.args.get("password")
    database = mod.args.get("database")
    db_url = make_mysql_url(
        host=host,
        port=port,
        user=user,
        password=password,
        database=database,
    )
    dbc = DBConnection(db_url)

    events = load_events(Path(mod.args.get("input_file")))
    with dbc.session_scope() as session:
        write_events(session, events)


if __name__ == "__main__":
    main()
