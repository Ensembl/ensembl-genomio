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

"""Download a Genbank file from NCBI from an accession.

Raises:
    DownloadError: if the download fails
"""


from os import PathLike
from pathlib import Path

import argschema
import requests


class DownloadError(Exception):
    """In case a download failed."""

    def __init__(self, msg):
        self.msg = msg


def download_genbank(accession: str, output_gb: PathLike) -> None:
    """
    Given a GenBank accession, download via NCBI Entrez service the corresponding file in GenBank format.

    Args:
        accession: INSDC Genbank record accession.
        output_gb: Path to the downloaded record in Genbank format.
    """

    # Get the list of assemblies for this accession
    entrez_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    entrez_params = {
        "db": "nuccore",
        "rettype": "gbwithparts",
        "retmode": "text",
    }
    entrez_params["id"] = accession
    result = requests.get(entrez_url, params=entrez_params, timeout=60)
    if result and result.status_code == 200:
        with Path(output_gb).open("wb") as gbff:
            gbff.write(result.content)
        print(f"GBF file write to {output_gb}")
        return
    raise DownloadError(f"Could not download the genbank ({accession}) file: {result}")


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    accession = argschema.fields.String(
        metadata={"required": True, "description": "Sequence accession required"}
    )
    output_file = argschema.fields.OutputFile(
        metadata={"required": True, "description": "Ouput Genbank path"}
    )


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    accession = mod.args["accession"]
    output_file = mod.args["output_file"]
    download_genbank(accession, output_file)


if __name__ == "__main__":
    main()
