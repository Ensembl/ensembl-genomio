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

"""Add or update one organism in an OSID server.

You can update an organism by setting a new pattern without setting the gene and transcript starts.
"""


import argparse
import json
from typing import Dict, Optional

import requests


def get_organism(url: str, user: str, key: str, org_data: Dict[str, str]) -> Optional[int]:
    """Retrieve the organism ID from OSID for a given organism.

    Returns None if no organism could be found with that name.
    """
    species = org_data["organismName"]

    body = {"organismName": species}

    get_org = url + "/organisms"
    result = requests.get(get_org, auth=(user, key), params=body, timeout=10)

    if result and result.status_code == 200:
        organisms = json.loads(result.content)

        if len(organisms) == 1:
            return organisms[0]["organismId"]
        if len(organisms) == 0:
            return None
        raise ValueError(f"Found several organisms for the name '{species}'")
    print(f"Error: {result}")
    return None


def add_organism(url: str, user: str, key: str, org_data: Dict[str, str]) -> None:
    """Add a new organism to OSID."""

    # Species is new, check we have all the keys
    expected = {"organismName", "template", "geneIntStart", "transcriptIntStart"}
    if set(expected) != set(org_data):
        raise ValueError(f"Missing values expected for a new organism: {org_data}")

    add_org_url = url + "/organisms"
    headers = {"Content-type": "application/json", "Accept": "application/json"}
    result = requests.post(add_org_url, auth=(user, key), headers=headers, json=org_data, timeout=10)

    if result and result.status_code == 200:
        print("Successfully added organism")
    else:
        print("Error: " + str(result))
        print(json.dumps(json.loads(result.content), indent=4))


def update_organism(url: str, user: str, key: str, org_data: Dict[str, str], organism_id: int) -> None:
    """Update an organism in OSID."""

    del org_data["organismName"]
    update_org_url = url + "/organisms" + "/" + str(organism_id)
    headers = {"Content-type": "application/json", "Accept": "application/json"}
    result = requests.put(update_org_url, auth=(user, key), headers=headers, json=org_data, timeout=10)

    if result and result.status_code == 204:
        print("Successfully updated organism")
    else:
        print("Error: " + str(result))
        print(json.dumps(json.loads(result.content), indent=4))


def add_or_update_organism(url: str, user, key: str, organism_data: Dict[str, str], update: bool) -> None:
    """Either add an organism or update it if it's already in OSID."""

    organism_id = get_organism(url, user, key, organism_data)
    if organism_id:
        print("Species already exist in the server")
        if update:
            print("Updating")
            update_organism(url, user, key, organism_data, organism_id)
        else:
            print("No changes: if you want to update it, add the --update argument")
    else:
        print("Adding species to the server")
        add_organism(url, user, key, organism_data)


def main():
    """Main script access point."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--url", type=str, required=True, help="OSID server url")
    parser.add_argument("--user", type=str, required=True, help="OSID user")
    parser.add_argument("--key", type=str, required=True, help="OSID authentification key")
    parser.add_argument(
        "--update", dest="update", action="store_true", help="Update the organism if it already exists"
    )

    parser.add_argument("--name", type=str, required=True, help="Species name to add")
    parser.add_argument("--template", type=str, required=True, help="Stable_id template")
    parser.add_argument("--gene_start", type=str, help="Start of gene numbering")
    parser.add_argument("--transcript_start", type=str, help="Start of transcript numbering")

    args = parser.parse_args()

    organism_data = {
        "organismName": args.name,
        "template": args.template,
    }
    if args.gene_start and args.transcript_start:
        organism_data["geneIntStart"] = args.gene_start
        organism_data["transcriptIntStart"] = args.transcript_start

    add_or_update_organism(args.url, args.user, args.key, organism_data, args.update)


if __name__ == "__main__":
    main()
