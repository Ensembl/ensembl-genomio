#!env python3
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



import argparse
from os import path
from typing import Any, List, Dict
from datetime import datetime

from ensembl.brc4.runnable.core_server import core_server


class UnsupportedEvent(Exception):
    pass


class StableIdEvent:
    def __init__(self, from_list: List[str] = [], to_list: List[str] = [], release: str = '', date: datetime = '') -> None:
        self.from_list = from_list
        self.to_list = to_list
        self.release = release
        self.date = date
        self.name = ""
        self._clean_lists()
    
    def __str__(self) -> str:
        from_str = ",".join(self.from_list)
        to_str = ",".join(self.to_list)
        release = self.release
        name = self.get_name()
        return f"From {from_str} to {to_str} = {name} in release {release}"
    
    def brc_format(self) -> List[str]:
        from_str = ",".join(self.from_list)
        to_str = ",".join(self.to_list)
        release = self.release
        date = self.date.strftime('%Y-%m')
        name = self.get_name()
        line_list = []
        for id in self.from_list:
            line = [
                id,
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
            new_id = self.to_list[0]
            line = (
                new_id,
                name,
                release,
                date,
                "",
                ""
            )
            line_list.append("\t".join(line))
        return line_list
    
    def _clean_lists(self):
        self.from_list = [id for id in self.from_list if id]
        self.to_list = [id for id in self.to_list if id]
    
    def add_from(self, from_id):
        self.from_list.append(from_id)
        self._clean_lists()
    
    def add_to(self, to_id):
        self.to_list.append(to_id)
        self._clean_lists()
    
    def add_release(self, release: str) -> None:
        self.release = release
    
    def add_date(self, date: datetime) -> None:
        self.date = date
    
    def name_event(self):
        if not self.from_list and len(self.to_list) == 1:
            self.name = "new"
        elif not self.to_list and len(self.from_list) == 1:
            self.name = "deletion"
        elif len(self.from_list) == 1 and len(self.to_list) == 1:
            self.name = "change"
        elif len(self.from_list) == 1 and len(self.to_list) > 1:
            self.name = "split"
        elif len(self.from_list) > 1 and len(self.to_list) == 1:
            self.name = "merge"
        elif len(self.from_list) > 1 and len(self.to_list) > 1:
            self.name = "mixed"
        else:
            raise UnsupportedEvent(f"Event {self.from_list} to {self.to_list} is not supported")
    
    def get_name(self):
        self.name_event()
        return self.name

class dump_stable_ids:

    def __init__(self, server: str) -> None:
        self.server = server
    
    def dump_history_json(self, output_file: str):
        pass

    def get_history(self) -> Any:
        
        sessions = self.get_mapping_sessions()

        events = []
        for session in sessions:
            print(f"Mapping session {session['release']}")
            pairs = self.get_pairs(session["id"])
            session_events = self.make_events(pairs)
            for event in session_events:
                event.add_release(session["release"])
                event.add_date(session["date"])
            events += session_events
        
        # Then analyse the pairs to make events
        return events
    
    def print_events(self, events: List[StableIdEvent], output_file: str) -> None:
        if not events:
            print("No events to print")
            return
        with open(output_file, 'w') as out_fh:
            for event in events:
                event_lines = event.brc_format()
                for line in event_lines:
                    out_fh.write(line + "\n")
    
    def get_mapping_sessions(self):
        query = """SELECT mapping_session_id, new_release, created
        FROM mapping_session
        """
        cursor = self.server.get_cursor()
        cursor.execute(query)

        sessions = []
        for db in cursor:
            date = db[2]
            session = {"id": db[0], "release": db[1], "date": date}
            sessions.append(session)
        return sessions

    def get_pairs(self, session_id):
        query = """SELECT old_stable_id, new_stable_id
        FROM stable_id_event
        WHERE (old_stable_id != new_stable_id OR old_stable_id IS NULL OR new_stable_id IS NULL)
            AND type="gene"
            AND mapping_session_id=%s
        """
        values = (session_id,)
        cursor = self.server.get_cursor()
        cursor.execute(query, values)

        pairs = []
        for db in cursor:
            pair = {"old_id": db[0], "new_id": db[1]}
            pairs.append(pair)
        print(f"{len(pairs)} stable id events")
        return pairs
    
    def make_events(self, pairs):

        from_list = {}
        to_list = {}
        for pair in pairs:
            old_id = pair["old_id"]
            new_id = pair["new_id"]

            if old_id in from_list:
                from_list[old_id].append(new_id)
            else:
                from_list[old_id] = [new_id]

            if new_id in to_list:
                to_list[new_id].append(old_id)
            else:
                to_list[new_id] = [old_id]
        
        # Remove empty elements
        for from_id in from_list:
            from_id_list = from_list[from_id]
            from_id_list = [id for id in from_id_list if id]
            from_list[from_id] = from_id_list
        for to_id in to_list:
            to_id_list = to_list[to_id]
            to_id_list = [id for id in to_id_list if id]
            to_list[to_id] = to_id_list
        
        events = []
        for old_id in from_list:
            if not old_id or old_id not in from_list: continue
            event = StableIdEvent([old_id], from_list[old_id])
            (event, from_list, to_list) = self.extend_event(event, from_list, to_list)
            events.append(event)
        
        # Remaining events should only be new genes
        for new_id in to_list:
            if not new_id: continue
            event = StableIdEvent(to_list[new_id], [new_id])
            events.append(event)
        
        stats = {}
        for event in events:
            name = event.get_name()
            if not name in stats:
                stats[name] = 1
            else:
                stats[name] += 1
        
        for stat in stats:
            print(f"\t{stat} = {stats[stat]}")

        return events
    
    def extend_event(self, event: StableIdEvent, from_list: Dict[str, List[str]], to_list: Dict[str, List[str]]):
        extended = True

        while(extended):
            extended = False
            
            # Extend the group in the to ids
            for to_id in event.to_list:
                if to_id in to_list:
                    to_from_ids = to_list[to_id]
                    # Add to the from list?
                    for to_from_id in to_from_ids:
                        if not to_from_id in event.from_list:
                            event.add_from(to_from_id)
                            extended = True
            
            # Extend the group in the from ids
            for from_id in event.from_list:
                if from_id in from_list:
                    from_to_ids = from_list[from_id]
                    # Add to the to list?
                    for from_to_id in from_to_ids:
                        if not from_to_id in event.to_list:
                            event.add_to(from_to_id)
                            extended = True
            
            # Clean up
            from_list = {from_id: from_list[from_id] for from_id in from_list if from_id not in event.from_list}
            to_list = {to_id: to_list[to_id] for to_id in to_list if to_id not in event.to_list}
        
        return (event, from_list, to_list)


def main():
    parser = argparse.ArgumentParser(description='Dump the stable ids history from a core db')

    parser.add_argument('--host', type=str, required=True, help='Server hostname')
    parser.add_argument('--port', type=str, required=True, help='Server port')
    parser.add_argument('--user', type=str, required=True, help='Server user')
    parser.add_argument('--password', type=str, help='Server password')
    parser.add_argument('--dbname', type=str, required=True, help='Database name')
    parser.add_argument('--output_file', type=str, required=True, help='Output file')

    # Optional
    args = parser.parse_args()

    # Start
    factory = core_server(host=args.host, port=args.port, user=args.user, password=args.password)
    factory.connect()
    factory.db.database = args.dbname
    dumper = dump_stable_ids(factory)
    events = dumper.get_history()
    dumper.print_events(events, args.output_file)


if __name__ == "__main__":
    main()