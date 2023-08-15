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

"""This script retrieves the organism metadata from an OSID server."""

import argparse
import sys
from typing import List
import json
import requests


class OSIDClient:
    """Client interface to an OSID server."""

    organisms_page = "organisms"

    def __init__(self, url: str, user: str, key: str) -> None:
        self.url = url
        self.user = user
        self.key = key

    def get_organism_data(self, species: str) -> List[dict]:
        """Retrieve the data for one or all organisms

        Args:
            species: optional, to only select one species instead of all of them

        Returns:
            List of dicts, one for each organism
        """
        organisms = []

        body = {
            "organismName": "*",
        }
        if species:
            body["organismName"] = species

        get_url = self.url + "/" + OSIDClient.organisms_page
        result = requests.get(get_url, auth=(self.user, self.key), params=body, timeout=30)

        if result and result.status_code == 200:
            organisms = json.loads(result.content)
        else:
            raise ValueError(f"Could not retrieve organism data from {get_url}")
        return organisms


def main():
    parser = argparse.ArgumentParser(description="Retrieve the list of organisms from an OSID server")
    parser.add_argument("--url", type=str, required=True, help="OSID server url")
    parser.add_argument("--user", type=str, required=True, help="OSID user")
    parser.add_argument("--key", type=str, required=True, help="OSID authentification key")
    parser.add_argument("--species", type=str, help="A specific species to check")
    args = parser.parse_args()

    # Get data
    client = OSIDClient(args.url, args.user, args.key)
    organisms = client.get_organism_data(args.species)

    # Display results
    print(f"{len(organisms)} Organisms in {args.url}:", file=sys.stderr)
    print(json.dumps(organisms, indent=4, sort_keys=True))


if __name__ == "__main__":
    main()
