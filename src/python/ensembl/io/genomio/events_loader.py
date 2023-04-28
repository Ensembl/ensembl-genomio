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
from typing import Dict, List
from dataclasses import dataclass

import argschema
from sqlalchemy.orm import Session

from ensembl.database import DBConnection
from ensembl.core.models import MappingSession, StableIdEvent


@dataclass
class IdEvent:
    """Simple representation for the events from the input file"""
    from_id: str
    to_id: str
    release: str
    release_date: str


class MapSession:
    """Simple mapping_sessions representation from the input file"""
    def __init__(self, release: str, release_date: str) -> None:
        self.release = release
        self.release_date = release_date
        self.events: List[IdEvent] = []

    def add_event(self, event: IdEvent) -> None:
        """Add an event to this mapping_session"""
        self.events.append(event)


def load_events(input_file: Path) -> List[IdEvent]:
    """Load events from input file.
    Expected tab file columns: old_id, new_id, event_name, release, release_date
    
    """
    events: List[IdEvent] = []

    with input_file.open("r") as events_fh:
        for line in events_fh:
            line.strip()
            if line == "":
                continue
            (from_id, to_id, _, release, release_date) = line.split("\t")
            event = IdEvent(from_id=from_id, to_id=to_id, release=release, release_date=release_date)
            events.append(event)
    return events


def write_events(session: Session, events: List[IdEvent], update: bool = False) -> None:
    """Insert the events in the core database.
    A mapping session is created for each different 'release'.
    
    """
    # First, create mapping_sessions based on the release
    mappings: Dict[str, MapSession] = {}
    for event in events:
        release = event.release

        if release in mappings:
            mappings[release].add_event(event)
        else:
            mappings[release] = MapSession(release, event.release_date)

    # Then, add the mapping, and the events for this mapping
    for release, mapping in mappings.items():
        print(f"Adding mapping for release {release} ({len(mapping.events)} events)")
        if update:
            map_session = MappingSession(new_release=mapping.release, created=mapping.release_date)
            session.add(map_session)
            session.flush()
            session.refresh(map_session)
            for event in mapping.events:
                id_event = StableIdEvent(
                    mapping_session_id=map_session.mapping_session_id,
                    old_stable_id=event.from_id,
                    new_stable_id=event.to_id,
                    id_type="gene",
                    old_version=1,
                    new_version=1,
                )
                session.add(id_event)
            session.commit()


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
    input_file = argschema.fields.InputFile(
        required=True,
        metadata={
            "description": (
                "Tab input events files in the format exported by the dumper:"
                "old_id, new_id, event_name, release, date"
            )
        },
    )
    update = argschema.fields.Boolean(
        default=False, required=False, metadata={"description": "Set this to actually make changes to the db"}
    )


def make_mysql_url(host: str, user: str, database: str, port: int = 0, password: str = "") -> str:
    """Create a mysql URL given server parameters for SQLalchemy."""
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
        write_events(session, events, mod.args["update"])


if __name__ == "__main__":
    main()
