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


from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Set, Tuple

import argschema

from ensembl.brc4.runnable.core_server import CoreServer


BRC4_START_DATE = datetime(2020, 5, 1)
IdsSet = Set[str]
DictToIdsSet = Dict[str, IdsSet]


class Pair:
    def __init__(self, old_id: Optional[str], new_id: Optional[str]) -> None:
        """Create a pair with an old_id and a new_id if provided"""

        if old_id is not None:
            self.old_id = old_id
        else:
            self.old_id = ""
        if new_id is not None:
            self.new_id = new_id
        else:
            self.new_id = ""

    def has_old_id(self) -> bool:
        return self.old_id != ""

    def has_new_id(self) -> bool:
        return self.new_id != ""

    def is_empty(self) -> bool:
        """Test if the current pair has no id."""

        return not (self.has_old_id() or self.has_new_id())


class UnsupportedEvent(ValueError):
    """If an event is not supported"""


class StableIdEvent:
    """Represents a stable id event from one gene set version to another one. Various events:
    - new genes
    - deleted genes
    - merged genes (several genes to one)
    - split genes (one gene to several)
    - mixed (several genes to several)

    Attributes:
        from_list: List of genes the previous gene set.
        to_list: List of genes in the new gene set.
        release: New gene set release name.
        date: Date of the new gene set.
        name: Name of the event (will be updated automatically).
        pairs: All pair of ids for this event.

    Any gene set before 2019-09 is dubbed pre-BRC4.

    """

    def __init__(
        self,
        from_list: Optional[Set[str]] = None,
        to_list: Optional[Set[str]] = None,
        release: Optional[str] = None,
        date: Optional[datetime] = None,
    ) -> None:
        if from_list is None:
            from_list = set()
        if to_list is None:
            to_list = set()
        self.from_set = self.clean_set(from_list)
        self.to_set = self.clean_set(to_list)
        self.release = release
        self.date = date
        self.name = ""
        self.pairs: List[Pair] = []

    def __str__(self) -> str:
        from_str = ",".join(self.from_set)
        to_str = ",".join(self.to_set)
        return f"From {from_str} to {to_str} = {self.get_name()} in release {self.release}"

    def brc_format_1(self) -> List[str]:
        """Returns a list events, one line per initial ID, in the following TSV format:
        - old gene id
        - event name
        - release
        - release date
        - list of old gene ids in the event (comma-separated)
        - list of new gene ids in the event (comma-separated)

        """
        from_str = ",".join(self.from_set)
        to_str = ",".join(self.to_set)
        release = self.get_full_release()
        if self.date:
            date = self.date.strftime("%Y-%m")
        else:
            date = "no_date"
        name = self.get_name()
        line_list = []
        for identifier in self.from_set:
            line = [
                identifier,
                name,
                release,
                date,
            ]
            if name in ("merge", "split", "mixed", "change"):
                line.append(from_str)
                line.append(to_str)
            else:
                line += ["", ""]
            line_list.append("\t".join(line))

        if self.get_name() == "new":
            new_id = [self.to_set][0]
            line = [new_id, name, release, date, "", ""]
            line_list.append("\t".join(line))
        return line_list

    def brc_format_2(self) -> List[str]:
        """Returns a list of combination of genes, one line per combination of old_id - new_ids, in the
        following TSV format:
        - old gene id
        - new gene id
        - event name
        - release
        - release date

        """
        release = self.get_full_release()
        if self.date:
            date = self.date.strftime("%Y-%m")
        else:
            date = "no_date"
        name = self.get_name()
        line_list = []

        for pair in self.pairs:
            line = [
                pair.old_id,
                pair.new_id,
                name,
                release,
                date,
            ]
            line_list.append("\t".join(line))
        return line_list

    @staticmethod
    def clean_set(this_list: Set) -> Set:
        """Removes any empty elements from a list.

        Args:
            this_list: list of items, so of which can be empty/None.

        Returns:
            The cleaned list.

        """
        return {identifier for identifier in this_list if identifier}

    def add_from(self, from_id: str) -> None:
        """Store an id in the from_set."""
        if from_id:
            self.from_set.add(from_id)

    def add_to(self, to_id: str) -> None:
        """Store an id in the from_set."""
        if to_id:
            self.to_set.add(to_id)

    def set_release(self, release: str) -> None:
        self.release = release

    def set_date(self, date: datetime) -> None:
        self.date = date

    def add_pair(self, pair: Pair) -> None:
        """Keeps a record of this pair.

        Args:
            pair: a Pair to record.

        Raises:
            ValueError: can't add an empty pair.

        """
        if pair.is_empty():
            raise ValueError(f"Expected at least one value in the given pair {pair}")
        self.pairs.append(pair)

    def get_full_release(self) -> str:
        """Returns the expanded release name, pre-BRC4 or `BRC4 = build`."""
        release = self.release
        date = self.date

        if date and date > BRC4_START_DATE:
            release = f"build {release}"
        else:
            release = f"pre-BRC4 {release}"

        return release

    def _name_event(self) -> None:
        """Identify the event name based on the old vs new id lists."""
        if not self.from_set and len(self.to_set) == 1:
            self.name = "new"
        elif not self.to_set and len(self.from_set) == 1:
            self.name = "deletion"
        elif len(self.from_set) == 1 and len(self.to_set) == 1:
            self.name = "change"
        elif len(self.from_set) == 1 and len(self.to_set) > 1:
            self.name = "split"
        elif len(self.from_set) > 1 and len(self.to_set) == 1:
            self.name = "merge"
        elif len(self.from_set) > 1 and len(self.to_set) > 1:
            self.name = "mixed"
        else:
            raise UnsupportedEvent(f"Event {self.from_set} to {self.to_set} is not supported")

    def clean_pairs(self) -> None:
        """Remove the empty old pairs when the event is not 'new'."""
        if not self.name:
            self._name_event()

        if self.name != "new":
            new_pairs = []
            for pair in self.pairs:
                if not pair.has_old_id():
                    continue
                new_pairs.append(pair)
            self.pairs = new_pairs

    def get_name(self) -> str:
        """Retrieve the name for this event, update it beforehand."""
        self._name_event()
        return self.name

    def add_pairs(self, pairs: List[Pair]) -> None:
        """Provided all the pairs, keep those that are used by this event.

        Args:
            pairs: list of Pair.

        """
        for pair in pairs:
            if (pair.has_old_id() and pair.old_id in self.from_set) or (
                pair.has_new_id() and pair.new_id in self.to_set
            ):
                # Core db contains an empty line to signify that an old id has been removed
                # in merge/split/mixed
                name = self.get_name()
                if (name != "deletion") and not pair.has_new_id():
                    continue
                self.add_pair(pair)


class DumpStableIDs:
    """An processor that create events from pairs of ids and can print those events out.

    Attributes:
        server: a core server set to a database, to retrieve the data from.

    """

    def __init__(self, server: CoreServer) -> None:
        self.server = server

    def get_history(self) -> List:
        """Retrieve all events from a database.

        Returns:
            A list of all events.

        """

        sessions = self.get_mapping_sessions()

        events = []
        for session in sessions:
            print(f"Mapping session {session['release']}")
            pairs = self.get_pairs(session["id"])
            session_events = self.make_events(pairs)
            for event in session_events:
                event.set_release(session["release"])
                event.set_date(session["date"])
            events += session_events

        # Then analyse the pairs to make events
        return events

    def print_events(self, events: List[StableIdEvent], output_file: Path) -> None:
        """Print events in a format for BRC.

        Args:
            events: list of events for a given genome.
            output_file: where the events will be printed.

        """
        if not events:
            print("No events to print")
            return
        with output_file.open("w") as out_fh:
            for event in events:
                event_lines = event.brc_format_2()
                for line in event_lines:
                    out_fh.write(line + "\n")

    def get_mapping_sessions(self) -> List[Dict]:
        """Retrieve the mapping sessions from the connected database.

        Returns:
            A list of sessions, as dicts: {'id: str, 'release': str, 'date': str}.

        """
        query = """SELECT mapping_session_id, new_release, created
        FROM mapping_session
        """
        cursor = self.server.get_cursor()
        cursor.execute(query)

        sessions: List[Dict[str, str]] = []
        for db in cursor:
            session = {"id": db[0], "release": db[1], "date": db[2]}
            sessions.append(session)
        return sessions

    def get_pairs(self, session_id: int) -> List[Pair]:
        """Retrieve all pair of ids for a given session.

        Args:
            session_id: id of a session from the connected database.

        Returns:
            A list of all pairs of ids, as dicts: {'old_id': str, 'new_id': str}.

        """
        query = """SELECT old_stable_id, new_stable_id
        FROM stable_id_event
        WHERE (old_stable_id != new_stable_id OR old_stable_id IS NULL OR new_stable_id IS NULL)
            AND type="gene"
            AND mapping_session_id=%s
        GROUP BY old_stable_id, new_stable_id, mapping_session_id
        """
        values = (session_id,)
        cursor = self.server.get_cursor()
        cursor.execute(query, values)

        pairs: List[Pair] = []
        for db in cursor:
            pair = Pair(old_id=db[0], new_id=db[1])
            pairs.append(pair)
        print(f"{len(pairs)} stable id events")
        return pairs

    def make_events(self, pairs: List[Pair]) -> List:
        """Given a list of pairs, create events.

        Args:
            pairs: list of Pair.

        Return:
            A list of events.

        """

        from_list: DictToIdsSet = {}
        to_list: DictToIdsSet = {}
        for pair in pairs:
            old_id = pair.old_id
            new_id = pair.new_id
            if old_id is None:
                old_id = ""
            if new_id is None:
                new_id = ""

            if old_id in from_list:
                from_list[old_id].add(new_id)
            else:
                from_list[old_id] = set([new_id])

            if new_id in to_list:
                to_list[new_id].add(old_id)
            else:
                to_list[new_id] = set([old_id])

        # Remove empty elements
        for from_id in from_list:
            from_list[from_id] = StableIdEvent.clean_set(from_list[from_id])
        for to_id in to_list:
            to_list[to_id] = StableIdEvent.clean_set(to_list[to_id])

        events: List[StableIdEvent] = []
        for old_id, from_old_list in from_list.items():
            if not old_id or old_id not in from_list:
                continue
            event = StableIdEvent(set([old_id]), set(from_old_list))
            (event, from_list, to_list) = self.extend_event(event, from_list, to_list)
            event.add_pairs(pairs)
            events.append(event)

        # Remaining events should only be new genes
        for new_id, to_new_list in to_list.items():
            if not new_id:
                continue
            event = StableIdEvent(set(to_new_list), set([new_id]))
            event.add_pairs(pairs)
            events.append(event)

        stats = {}
        for event in events:
            name = event.get_name()
            event.clean_pairs()
            if name not in stats:
                stats[name] = 1
            else:
                stats[name] += 1

        for stat, value in stats.items():
            print(f"\t{stat} = {value}")

        return events

    def extend_event(
        self, event: StableIdEvent, from_list: DictToIdsSet, to_list: DictToIdsSet
    ) -> Tuple[StableIdEvent, DictToIdsSet, DictToIdsSet]:
        """Given an event, aggregate ids in the 'from' and 'to' sets, to connect the whole group.

        Args:
            event: the event to extend.
            from_list: A dict a the from ids, and their corresponding to ids.
            to_list: A dict of the to ids, and their corresponding from ids.

        Returns:
            A tuple of the extended event, and the from_list and to_list from which the ids that
            have been added to the event have been removed.

        """

        extended = True

        while extended:
            extended = False

            # Extend the group in the to ids
            for to_id in event.to_set:
                if to_id in to_list:
                    to_from_ids: IdsSet = to_list[to_id]
                    # Add to the from list?
                    for to_from_id in to_from_ids:
                        if to_from_id not in event.from_set:
                            event.add_from(to_from_id)
                            extended = True

            # Extend the group in the from ids
            for from_id in event.from_set:
                if from_id in from_list:
                    from_to_ids = from_list[from_id]
                    # Add to the to list?
                    for from_to_id in from_to_ids:
                        if from_to_id not in event.to_set:
                            event.add_to(from_to_id)
                            extended = True

        # Clean up
        from_list = {from_id: from_list[from_id] for from_id in from_list if from_id not in event.from_set}
        to_list = {to_id: to_list[to_id] for to_id in to_list if to_id not in event.to_set}

        return (event, from_list, to_list)


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
    output_file = argschema.fields.OutputFile(required=True, metadata={"description": "Output file"})


def main() -> None:
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    args = mod.args

    # Start
    factory = CoreServer(
        host=args.get("host"), port=args.get("port"), user=args.get("user"), password=args.get("password")
    )
    factory.set_database(args.get("database"))
    dumper = DumpStableIDs(factory)
    events = dumper.get_history()
    dumper.print_events(events, Path(args.get("output_file")))


if __name__ == "__main__":
    main()
