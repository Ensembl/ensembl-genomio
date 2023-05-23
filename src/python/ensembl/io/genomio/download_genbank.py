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


import requests
import argschema


class DownloadError(Exception):
    """In case a download failed."""

    def __init__(self, msg):
        self.msg = msg


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    accession = argschema.fields.String(
        metadata={"required": True, "description": "Sequence accession required"}
    )


def download_genbank(accession: str) -> str:
    """
    Given a GenBank accession, download the corresponding file in GenBank format
    """
    dl_file = f"{accession}.gb"

    # Get the list of assemblies for this accession
    e_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    e_params = {
        "db": "nuccore",
        "rettype": "gbwithparts",
        "retmode": "text",
    }
    e_params["id"] = accession
    result = requests.get(e_url, params=e_params, timeout=60)
    if result and result.status_code == 200:
        with open(dl_file, "wb") as gff:
            gff.write(result.content)
        print(f"GFF file write to {dl_file}")
        return dl_file
    raise DownloadError(f"Could not download the genbank file: {result}")


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    accession = mod.args["accession"]
    download_genbank(accession)


if __name__ == "__main__":
    main()
