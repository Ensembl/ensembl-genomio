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
"""Obtain and record the assembly status for a set of INSDC accession(s) using NCBI's datasets CLI tool."""

__all__ = [
    "extract_assembly_metadata",
    "fetch_datasets_reports",
    "fetch_accessions_from_core_dbs",
    "generate_report_tsv",
    "get_assembly_accessions",
    "singularity_image_setter",
]

import csv
from dataclasses import dataclass
import json
import logging
import os
from pathlib import Path
import re

from spython.main import Client
from sqlalchemy.engine import URL
from sqlalchemy import text

import ensembl.io.genomio
from ensembl.io.genomio.utils.json_utils import print_json
from ensembl.utils import StrPath
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.database import DBConnection
from ensembl.utils.logging import init_logging_with_args


DATASETS_SINGULARITY = {
    "datasets_version_url": "docker://ensemblorg/datasets-cli:latest",
}


class UnsupportedFormatError(Exception):
    """When a string does not have the expected format."""


@dataclass
class ReportStructure:
    """Stores key report meta information."""

    species_name: str = ""
    taxon_id: int = 0
    strain: str = "NA"
    assembly_name: str = ""
    assembly_type: str = ""
    accession: str = ""
    paired_assembly: str = "NA"
    last_updated: str = ""
    assembly_status: str = "NA"
    assembly_notes: str = "NA"

    def to_dict(self) -> dict[str, str]:
        """Returns a dictionary representation of this object."""
        return {
            "Species Name": self.species_name,
            "Taxon ID": str(self.taxon_id),
            "Isolate/Strain": self.strain,
            "Asm name": self.assembly_name,
            "Assembly type": self.assembly_type,
            "Asm accession": self.accession,
            "Paired assembly": self.paired_assembly,
            "Asm last updated": self.last_updated,
            "Asm status": self.assembly_status,
            "Asm notes": self.assembly_notes,
        }

    def header(self) -> list[str]:
        """Returns the dictionary keys matching each of the properties of the report."""
        return list(self.to_dict().keys())

    def values(self) -> list[str]:
        """Returns the values of each of the properties of the report."""
        return list(self.to_dict().values())


def singularity_image_setter(sif_cache_dir: Path | None, datasets_version: str | None) -> Client:
    """Parse ENV and User specified variables related to `datasets` singularity SIF
    container and define version and location of container.

    Args:
        sif_cache_dir: Path to locate existing, or download new SIF container image.
        datasets_version: URL of singularity container (custom `datasets` version if desired).

    Returns:
        `spython.main.client` instance of singularity container image housing `datasets`.
    """

    # Set singularity cache dir from user defined path or use environment
    if sif_cache_dir and sif_cache_dir.is_dir():
        image_dl_path = sif_cache_dir
        logging.info(f"Using user-defined cache_dir: '{image_dl_path}'")
    elif os.environ.get("NXF_SINGULARITY_CACHEDIR"):
        image_dl_path = Path(os.environ["NXF_SINGULARITY_CACHEDIR"])
        logging.info(
            f"Using preferred nextflow singularity cache dir 'NXF_SINGULARITY_CACHEDIR': {image_dl_path}"
        )
    elif os.environ.get("SINGULARITY_CACHEDIR"):
        image_dl_path = Path(os.environ["SINGULARITY_CACHEDIR"])
        logging.info(
            f"Using the default singularity installation cache dir 'SINGULARITY_CACHEDIR': {image_dl_path}"
        )
    else:
        image_dl_path = Path()
        logging.warning(f"Unable to set singularity cache dir properly, using CWD {image_dl_path}")

    # Set the datasets version URL
    if datasets_version is None:
        container_url = DATASETS_SINGULARITY["datasets_version_url"]
        logging.info(f"Using default 'ncbi datasets' version '{container_url}'")
    else:
        container_url = datasets_version
        logging.info(f"Using user defined 'ncbi datasets' version '{container_url}'")

    # Pull or load pre-existing 'datasets' singularity container image.
    datasets_image = Client.pull(container_url, stream=False, pull_folder=image_dl_path, quiet=True)

    return datasets_image


def get_assembly_accessions(src_file: StrPath) -> list[str]:
    """Returns the list of assembly accessions found in the provided file.

    Args:
        src_file: Path to file with one line per INSDC assembly accession.

    Raises:
        UnsupportedFormatError: If an accession does not match the INSDC assembly accession format.
    """
    query_accessions: list[str] = []
    with Path(src_file).open(mode="r") as fin:
        for line in fin.readlines():
            line = line.strip()
            match = re.match(r"^GC[AF]_[0-9]{9}\.[1-9][0-9]*$", line)
            if not match:
                raise UnsupportedFormatError(f"Could not recognize GCA/GCF accession format: {line}")
            query_accessions.append(line)
    return query_accessions


def fetch_accessions_from_core_dbs(src_file: StrPath, server_url: URL) -> dict[str, str]:
    """Obtain the associated INSDC accession given a set of core database names and a database server URL.

    The accession information is obtained from the `meta` table's meta key `assembly.accession`.

    Args:
        src_file: File path with list of core database names.
        server_url: Database server URL.

    Returns:
        Dict of core database names (key) and their corresponding INSDC assembly accession (value).
    """

    core_accn_meta = {}
    database_count = 0
    count_accn_found = 0

    with Path(src_file).open("r") as fin:
        for line in fin.readlines():
            core_db = line.strip()
            database_count += 1
            db_connection_url = server_url.set(database=core_db)
            db_connection = DBConnection(db_connection_url)
            with db_connection.begin() as conn:
                query_result = conn.execute(
                    text('SELECT meta_value FROM meta WHERE meta_key = "assembly.accession";')
                ).fetchall()

            if not query_result:
                logging.warning(f"No accessions found in core: {core_db}")
            elif len(query_result) == 1:
                count_accn_found += 1
                asm_accession = query_result.pop()[0]
                logging.info(f"{core_db} -> assembly.accession[{asm_accession}]")
                core_accn_meta[core_db] = asm_accession
            else:
                logging.warning(f"Core {core_db} has {len(query_result)} assembly.accessions")

    logging.info(
        f"From initial input core databases ({database_count}), obtained ({count_accn_found}) accessions"
    )

    return core_accn_meta


def fetch_datasets_reports(
    sif_image: Client, assembly_accessions: dict[str, str], download_directory: StrPath, batch_size: int
) -> dict[str, dict]:
    """Obtain assembly reports in JSON format for each assembly accession via `datasets` CLI.

    Args:
        sif_image: Instance of `Client.loaded()` singularity image.
        assembly_accessions: Dictionary of accession source <> assembly accessions pairs.
        download_directory: Directory path to store assembly report JSON files.
        batch_size: Number of assembly accessions to batch submit to `datasets`.

    Returns:
        Dictionary of accession source and its associated assembly report.

    Raises:
        ValueError: If result returned by `datasets` is not a string.
        RuntimeError: If there was an error raised by `datasets`.

    """
    master_accn_list = list(assembly_accessions.values())
    combined_asm_reports = {}

    # Setting the number of combined accessions to query in a single call to datasets
    list_split = list(range(0, len(master_accn_list), batch_size))
    accn_subsample = [master_accn_list[ind : ind + batch_size] for ind in list_split]

    datasets_command = ["datasets", "summary", "genome", "accession"]
    for accessions in accn_subsample:
        # Make call to singularity datasets providing a multi-accession query
        client_return = Client.execute(
            image=sif_image, command=datasets_command + accessions, return_result=True, quiet=True
        )
        raw_result = client_return["message"]

        ## Test what result we have obtained following execution of sif image and accession value
        # Returned a list, i.e. datasets returned a result to client.execute
        # Returned a str, i.e. no datasets result obtained exited with fatal error
        if isinstance(raw_result, list):
            result = raw_result[0]
        else:
            result = raw_result
        if not isinstance(result, str):
            raise ValueError("Result obtained from datasets is not a string")
        if re.search("^FATAL", result):
            raise RuntimeError(f"Singularity image execution failed! -> '{result.strip()}'")

        tmp_asm_dict = json.loads(result)
        if not tmp_asm_dict["total_count"]:
            logging.warning(f"No assembly report found for accession(s) {accessions}")
            continue

        logging.info(f"Assembly report obtained for accession(s) {accessions}")
        batch_reports_json = tmp_asm_dict["reports"]
        for assembly_report in batch_reports_json:
            accession = assembly_report["accession"]
            asm_json_outfile = Path(download_directory, f"{accession}.asm_report.json")
            print_json(asm_json_outfile, assembly_report)
            # Save assembly report into source key<>report dict
            for src_key, accession_core in assembly_accessions.items():
                if accession == accession_core:
                    combined_asm_reports[src_key] = assembly_report

    return combined_asm_reports


def extract_assembly_metadata(assembly_reports: dict[str, dict]) -> dict[str, ReportStructure]:
    """Parse assembly reports and extract specific key information on status and related fields.

    Args:
        assembly_reports: Key value pair of source name <> assembly report.

    Returns:
        Parsed assembly report meta (source, meta).
    """
    parsed_meta = {}

    for source, asm_report in assembly_reports.items():
        asm_meta_info = ReportStructure()

        # Mandatory meta key parsing:
        asm_meta_info.accession = asm_report["accession"]
        asm_meta_info.assembly_name = asm_report["assembly_info"]["assembly_name"]
        asm_meta_info.assembly_type = asm_report["assembly_info"]["assembly_type"]
        asm_meta_info.assembly_status = asm_report["assembly_info"]["assembly_status"]
        asm_meta_info.species_name = asm_report["organism"]["organism_name"]
        asm_meta_info.taxon_id = int(asm_report["organism"]["tax_id"])

        ## Non-mandatory meta key parsing:
        assembly_meta_keys = asm_report["assembly_info"].keys()
        organism_keys = asm_report["organism"].keys()

        # check for genome_notes:
        if "genome_notes" in assembly_meta_keys:
            complete_notes = ", ".join(asm_report["assembly_info"]["genome_notes"])
            asm_meta_info.assembly_notes = complete_notes

        # check for biosample:
        if "biosample" in assembly_meta_keys:
            asm_meta_info.last_updated = asm_report["assembly_info"]["biosample"]["last_updated"]

        # check for paired assembly:
        if "paired_assembly" in assembly_meta_keys:
            asm_meta_info.paired_assembly = asm_report["assembly_info"]["paired_assembly"]["accession"]

        # check for isolate/strain type:
        if "infraspecific_names" in organism_keys:
            organism_type_keys = asm_report["organism"]["infraspecific_names"].keys()
            if "isolate" in organism_type_keys:
                asm_meta_info.strain = asm_report["organism"]["infraspecific_names"]["isolate"]
            elif "strain" in organism_type_keys:
                asm_meta_info.strain = asm_report["organism"]["infraspecific_names"]["strain"]

        parsed_meta[source] = asm_meta_info

    return parsed_meta


def generate_report_tsv(
    parsed_asm_reports: dict[str, ReportStructure],
    query_type: str,
    output_directory: StrPath = Path(),
    outfile_name: str = "AssemblyStatusReport",
) -> None:
    """Generate and write the assembly report to a TSV file.

    Args:
        parsed_asm_reports: Parsed assembly report meta.
        query_type: Type of query (either core databases or accessions).
        output_directory: Directory to store report TSV file.
        outfile_name: Name to give to the output TSV file.
    """
    tsv_outfile = Path(output_directory, f"{outfile_name}.tsv")

    header_list = next(iter(parsed_asm_reports.values())).header()
    header_list = [query_type.capitalize().replace("_", " ")] + header_list

    with open(tsv_outfile, "w+") as tsv_out:
        writer = csv.writer(tsv_out, delimiter="\t", lineterminator="\n")
        writer.writerow(header_list)
        for core, report_meta in parsed_asm_reports.items():
            final_asm_report = [core] + report_meta.values()
            writer.writerow(final_asm_report)


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    # Create parser with common arguments to be used by both subparsers
    base_parser = ArgumentParser(add_help=False)
    base_parser.add_argument_dst_path(
        "--reports_dir",
        default=Path("assembly_report_jsons"),
        help="path to folder where the assembly report JSON files is stored",
    )
    base_parser.add_argument(
        "--assembly_report_name",
        metavar="NAME",
        default="AssemblyStatusReport",
        help="file name used for the assembly report TSV output file",
    )
    base_parser.add_argument(
        "--datasets_version_url",
        type=str,
        metavar="URL",
        help="datasets version, e.g. docker://ensemblorg/datasets-cli:latest",
    )
    base_parser.add_argument_src_path(
        "--cache_dir",
        default=Path(os.environ.get("NXF_SINGULARITY_CACHEDIR", "")),
        metavar="SINGULARITY_CACHE",
        help="folder path to user generated singularity container housing NCBI tool 'datasets'",
    )
    base_parser.add_numeric_argument(
        "--datasets_batch_size",
        type=int,
        min_value=1,
        default=100,
        metavar="BATCH_SIZE",
        help="number of accessions requested in one query to datasets",
    )
    base_parser.add_log_arguments(add_log_file=True)
    # Add subparsers with their parent being the base parser with the common arguments
    subparsers = parser.add_subparsers(title="report assembly status from", required=True, dest="src")
    # Specific arguments required when using Ensembl core database names as source
    core_db_parser = subparsers.add_parser(
        "core_db", parents=[base_parser], help="list of Ensembl core databases"
    )
    core_db_parser.add_argument_src_path(
        "--input",
        required=True,
        help="file path with list of Ensembl core database(s) to retrieve query accessions from",
    )
    core_db_parser.add_server_arguments()
    # Specific arguments required when using assembly accessions as source
    accessions_parser = subparsers.add_parser(
        "accession", parents=[base_parser], help="list of INSDC accessions"
    )
    accessions_parser.add_argument_src_path(
        "--input", required=True, help="file path with list of assembly INSDC query accessions"
    )

    args = parser.parse_args()
    init_logging_with_args(args)

    # Get accessions on cores list or use user accession list directly
    if args.src == "core_db":
        query_accessions = fetch_accessions_from_core_dbs(args.input, args.url)
    else:
        query_accessions = {x: x for x in get_assembly_accessions(args.input)}

    # Parse singularity setting and define the SIF image for 'datasets'
    datasets_image = singularity_image_setter(args.cache_dir, args.datasets_version_url)

    # Datasets query implementation for one or more batched accessions
    assembly_reports = fetch_datasets_reports(
        datasets_image, query_accessions, args.reports_dir, args.datasets_batch_size
    )

    # Extract the key assembly report meta information for reporting status
    key_assembly_report_meta = extract_assembly_metadata(assembly_reports)

    # Produce the finalized assembly status report TSV from set of parsed 'datasets summary report'
    generate_report_tsv(key_assembly_report_meta, args.src, args.reports_dir, args.assembly_report_name)
