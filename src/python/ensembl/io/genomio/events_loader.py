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

"""Provided a file with events, load them in a core database.

cf the load_events functions for the events tab file format.
"""


from dataclasses import dataclass
from os import PathLike
from pathlib import Path
import re
from typing import Dict, Generator, List, Tuple

import argschema
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from ensembl.database import DBConnection
from ensembl.core.models import MappingSession, StableIdEvent


@dataclass
class IdEvent:
    """Simple representation for the events from the input file"""

    from_id: str
    to_id: str
    event: str
    release: str
    release_date: str

    def __str__(self) -> str:
        fields = [self.from_id, self.to_id, self.event, self.release, self.release_date]
        return "\t".join(fields)


class MapSession:
    """Simple mapping_sessions representation from the input file"""

    def __init__(self, release: str, release_date: str) -> None:
        self.release = release
        self.release_date = release_date
        self.events: List[IdEvent] = []

    def add_event(self, event: IdEvent) -> None:
        """Add an event to this mapping_session"""
        self.events.append(event)


class EventCollection:
    """ "Collection of events with loader/writer in various formats."""

    def __init__(self) -> None:
        self.events: List[IdEvent] = []

    def load_events(self, input_file: PathLike):
        """Load events from input file.
        Expected tab file columns: old_id, new_id, event_name, release, release_date

        """
        events: List[IdEvent] = []

        with Path(input_file).open("r") as events_fh:
            for line in events_fh:
                line.strip()
                if line == "":
                    continue
                (from_id, to_id, event_name, release, release_date) = line.split("\t")
                event = IdEvent(
                    from_id=from_id, to_id=to_id, event=event_name, release=release, release_date=release_date
                )
                events.append(event)
        self.events = events

    def add_deletes(
        self, genes: List[str], release_name: str = "release_name", release_date: str = "release_date"
    ) -> None:
        """Add deletion events from a list of deleted genes."""
        for gene_id in genes:
            event = IdEvent(
                from_id=gene_id, to_id="", event="deletion", release=release_name, release_date=release_date
            )
            self.events.append(event)

    def load_events_from_gene_diff(
        self, input_file: PathLike, release_name: str = "release_name", release_date: str = "release_date"
    ):
        """Load events from input file from gene_diff."""
        event_name = "event"
        loaded_event = set()

        with Path(input_file).open("r") as events_fh:
            for line in events_fh:
                if line.startswith("//") or line == "":
                    continue
                (_, event_string, _) = line.split("\t")
                for pair in self._parse_gene_diff_event(event_string):
                    (from_id, to_id) = pair
                    fingerprint = f"{from_id} {to_id}"
                    if fingerprint in loaded_event:
                        print(f"Duplicated event, skipped: {fingerprint}")
                        continue
                    loaded_event.add(fingerprint)
                    event = IdEvent(
                        from_id=from_id,
                        to_id=to_id,
                        event=event_name,
                        release=release_name,
                        release_date=release_date,
                    )
                    self.events.append(event)

    def _parse_gene_diff_event(self, event_string: str) -> Generator[Tuple[str, str], None, None]:
        """Gets all the pairs of IDs from an event string from gene diff."""
        splitter = re.compile(r"(~|\+|=\+|=-|=!|=|>|<)")
        parts = re.split(splitter, event_string)
        if len(parts) != 3:
            return
        [from_ids, sep, to_ids] = parts

        # Identical gene: no need to keep in the history
        if sep == "~":
            return

        for from_id in from_ids.split(":"):
            for to_id in to_ids.split(":"):
                yield (from_id, to_id)

    def remap_to_ids(self, map_dict: Dict[str, str]):
        """Using a mapping dict, remap the to_id of all events."""

        no_map = 0
        for event in self.events:
            if not event.to_id:
                continue
            if event.to_id in map_dict:
                event.to_id = map_dict[event.to_id]
            else:
                print(f"No map for to_id {event.to_id}")
                no_map += 1

        if no_map:
            raise ValueError(f"No map for {no_map} event to_ids")

    def write_events_to_file(self, output_file: PathLike) -> None:
        """Write the events to a file."""
        with Path(output_file).open("w") as out_fh:
            print(f"Write {len(self.events)} events to {output_file}")
            for event in self.events:
                out_fh.write(f"{event}\n")

    def write_events_to_db(self, session: Session, update: bool = False) -> None:
        """Insert the events in the core database.
        A mapping session is created for each different 'release'.

        """
        # First, create mapping_sessions based on the release
        mappings: Dict[str, MapSession] = {}
        for event in self.events:
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


def main() -> None:
    """Main entrypoint"""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)

    # Start
    db_url = URL.create(
        "mysql",
        mod.args["user"],
        mod.args.get("password"),
        mod.args["host"],
        mod.args["port"],
        mod.args.get("database"),
    )
    dbc = DBConnection(db_url)

    collection = EventCollection()
    collection.load_events(mod.args.get("input_file"))

    with dbc.session_scope() as session:
        collection.write_events_to_db(session, mod.args["update"])


if __name__ == "__main__":
    main()
