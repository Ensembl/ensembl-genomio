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
"""Download a Genbank file from NCBI from an accession."""

__all__ = ["DownloadError", "download_genbank"]

import logging
from os import PathLike
from pathlib import Path

import requests

import ensembl.io.genomio
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


class DownloadError(Exception):
    """In case a download failed."""

    def __init__(self, msg: str) -> None:
        self.msg = msg


def download_genbank(accession: str, output_file: PathLike) -> None:
    """Given a GenBank accession, download the corresponding file in GenBank format.

    Uses NCBI Entrez service to fetch the data.

    Args:
        accession: INSDC Genbank record accession.
        output_file: Path to the downloaded record in Genbank format.

    Raises:
        DownloadError: If the download fails.

    """

    # Get the list of assemblies for this accession
    entrez_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    entrez_params = {
        "db": "nuccore",
        "rettype": "gbwithparts",
        "retmode": "text",
    }
    entrez_params["id"] = accession
    logging.debug(f"Getting file from {entrez_url} with params {entrez_params}")
    result = requests.get(entrez_url, params=entrez_params, timeout=60)
    if result and result.status_code == 200:
        with Path(output_file).open("wb") as gbff:
            gbff.write(result.content)
        logging.info(f"GenBank file written to {output_file}")
        return
    raise DownloadError(f"Could not download the genbank ({accession}) file: {result}")


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description="Download a sequence from GenBank.")
    parser.add_argument("--accession", required=True, help="Sequence accession")
    parser.add_argument_dst_path("--output_file", required=True, help="Output GenBank file")
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    download_genbank(accession=args.accession, output_file=args.output_file)


if __name__ == "__main__":
    main()
