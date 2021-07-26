#!env python3

import argparse
import requests
import json

def get_organisms(url, user, key, species):

    body = {
            "organismName": "*",
            }
    if species:
        body["organismName"] = species
    
    get_org = url + "/organisms"
    result = requests.get(get_org, auth=(user, key), params=body)

    if result and result.status_code == 200:
        organisms = json.loads(result.content)
        count = len(organisms)
        print("%d Organisms in %s:" % (count, get_org))
        print(json.dumps(organisms, indent=4, sort_keys=True))
    else:
        print("Error: " + str(result))


def main():
    parser = argparse.ArgumentParser(description='Retrieve the list of organisms from an OSID server')
    
    parser.add_argument('--url', type=str, required=True,
                help='OSID server url')
    parser.add_argument('--user', type=str, required=True,
                help='OSID user')
    parser.add_argument('--key', type=str, required=True,
                help='OSID authentification key')
    parser.add_argument('--species', type=str,
                help='A specific species to check')

    args = parser.parse_args()
    
    get_organisms(args.url, args.user, args.key, args.species)

if __name__ == "__main__":
    main()
