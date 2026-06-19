# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Update the switcher JSON file with a new version entry.

Adds a new version to the PyData Sphinx Theme version switcher dropdown.
Usage: `python update_switcher.py --version <version> --switcher <path_to_switcher_json> \
            [--base-url <base_url>]`
"""

import json

from ensembl.utils import argparse


def main() -> None:
    """Script's main entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True, help="Version to add to the switcher")
    parser.add_argument(
        "--base-url",
        default="https://ensembl.github.io/ensembl-genomio",
        help="Base URL of the documentation site",
    )
    parser.add_argument_src_path("--switcher", required=True, help="Path to the switcher JSON file")
    args = parser.parse_args()
    # Avoid duplicates and remove "(latest)" from last version
    versions = json.loads(args.switcher.read_text()) if args.switcher.exists() else []
    versions = [v for v in versions if v["version"] != args.version]
    if versions:
        versions[0]["name"] = versions[0]["version"]
    new_entry = {
        "name": f"{args.version} (latest)",
        "version": args.version,
        "url": f"{args.base_url.rstrip('/')}/{args.version}/",
    }
    versions.insert(0, new_entry)
    args.switcher.write_text(json.dumps(versions, indent=2))


if __name__ == "__main__":
    main()
