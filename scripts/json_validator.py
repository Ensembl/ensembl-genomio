#!env python3

import argparse
import json
from jsonschema import validate


def run(json_file, json_schema):

    file_data = get_json(json_file)
    schema_data = get_json(json_schema)

    validate(instance=file_data, schema=schema_data)

def get_json(json_path):
    with open(json_path) as json_file:
        data = json.load(json_file)
        return data


def main():
    parser = argparse.ArgumentParser(description='Validate a json file against a json schema')
    parser.add_argument('json_file', metavar='json_file', type=str,
                help='the json file to test')
    parser.add_argument('json_schema', metavar='json_schema', type=str,
                help='the json schema to test against')

    args = parser.parse_args()
    run(args.json_file, args.json_schema)

if __name__ == "__main__":
    main()
