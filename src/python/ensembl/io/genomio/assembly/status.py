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
    "check_parameterization",
    "datasets_asm_reports",
    "extract_assembly_metadata",
    "fetch_accessions_from_cores",
    "generate_report_tsv",
    "resolve_query_type",
    "singularity_image_setter",
]

import csv
import json
import logging
import os
from os import PathLike
from pathlib import Path
import re
from typing import Dict, List, Tuple, Union

from spython.main import Client
from sqlalchemy.engine import URL
from sqlalchemy import text

from ensembl.io.genomio.utils.json_utils import print_json
from ensembl.io.genomio.database.dbconnection_lite import DBConnectionLite as dbc
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


DATASETS_SINGULARITY = {
    "datasets_version_url": "docker://ensemblorg/datasets-cli:latest",
}


class UnsupportedFormatError(Exception):
    """When a string does not have the expected format."""


class ReportStructure(dict):
    """Dict setter class of key report meta information"""

    def __init__(self):
        dict.__init__(self)
        self.update(
            {
                "Species Name": "",
                "Taxon ID": "",
                "Isolate/Strain": "",
                "Asm name": "",
                "Assembly type": "",
                "Asm accession": "",
                "Paired assembly": "",
                "Asm last updated": "",
                "Asm status": "",
                "Asm notes": "",
            }
        )


def singularity_image_setter(sif_cache_dir: Path, datasets_version: str) -> Client:
    """Parse ENV and User specified variables related to 'datasets' singularity SIF
    container and define version and location of container.

    Args:
        sif_cache_dir: Path to locate existing, or download new SIF container image.
        datasets_version: URL of singularity container (custom 'datasets' version if desired)

    Returns:
        `spython.main.client` instance of singularity container image housing 'datasets'.
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


def check_parameterization(input_cores: Path, input_accessions: Path, db_host: str, db_port: int) -> Path:
    """Detect the kind of user input (cores/accessions) and determine any missing or
    incorrect parameterization.

    Args:
        input_cores: Input core(s) list file name.
        input_accessions: Input accession (s) list file name.
        db_host: Host server name
        db_port: Host server port

    Returns:
        User input file used in assembly status querying
    """
    # Input core names centered run
    if input_cores:
        logging.info(f"Performing assembly status report using core db list file: {input_cores}")
        if db_host is None or db_port is None:
            raise RuntimeError("Core database names require both arguments '--host' and '--port'")
        return input_cores
    # Accession centered run
    logging.info(f"Performing assembly status report using INSDC accession list file: {input_accessions}")
    return input_accessions


def resolve_query_type(
    query_list: list, partial_url: URL, input_cores: str, input_accessions: str
) -> Union[Tuple[Dict, str]]:
    """Function to identify the kind of queries being passed by user,
    then extract the queries (core names or accessions) and store each with appropriate identifier.

    Args:
        query_list: List of user defined queries either core names, or accessions
        partial_url: A partial MYSQL connection URL (host:port)
        input_cores: Arg parse param '--input_cores'
        input_accessions: Arg parse param '--input_accessions'

    Returns:
        User queries stored as identifier[(core db name | UniqueID#)] : accession
    """

    query_accessions: Dict = {}
    query_type: str = ""

    if input_cores and input_accessions is None:
        query_accessions = fetch_accessions_from_cores(query_list, partial_url)
        query_type = "CoreDB"
    elif input_cores is None and input_accessions:
        query_type = "Accession"
        for accession in query_list:
            match = re.match(r"(GC[AF])_([0-9]{3})([0-9]{3})([0-9]{3})\.?([0-9]+)", accession)
            if not match:
                raise UnsupportedFormatError(f"Could not recognize GCA accession format: {accession}")
            query_accessions[accession] = accession

    return query_accessions, query_type


def fetch_accessions_from_cores(database_names: List, connection_url: URL) -> Dict:
    """Obtain the associated INSDC accession [meta.assembly.accession] given a set of core(s) names
    and a MYSQL server host.

    Args:
        database_names: Set of names for one or more core databases
        connection_url: Partial MYSQL host name : port

    Returns:
        Dict of core name(s) (key) and its INSDC assembly.accession (value)
    """

    core_accn_meta = {}
    core_list_count = len(database_names)
    count_accn_found = 0

    for core in database_names:
        db_connection_url = connection_url.set(database=core)
        db_connection = dbc(f"{db_connection_url}")
        with db_connection.connect() as conn:
            query_result = conn.execute(
                text('SELECT meta_value FROM meta WHERE meta_key = "assembly.accession";')
            ).fetchall()

        if query_result is None:
            logging.warning(f"No accessions found in core: {core}")
        elif len(query_result) == 1:
            count_accn_found += 1
            asm_accession = query_result.pop()[0]
            logging.info(f"{core} -> assembly.accession[{asm_accession}]")
            core_accn_meta[core] = asm_accession
        else:
            logging.warning(f"Core {core} has {len(query_result)} assembly.accessions")

    logging.info(f"From initial input cores ({core_list_count}), obtained ({count_accn_found}) accessions")

    return core_accn_meta


def datasets_asm_reports(
    sif_image: str, assembly_accessions: dict, download_directory: PathLike, batch_size: int
) -> Dict:
    """Obtain assembly report(s) JSONs for one or more queries made to datasets CLI.

    Args:
        sif_image: Instance of Client.loaded singularity image.
        assembly_accessions: Dict of core accessions.
        download_directory: Dir path to store assembly report JSON files.
        batch_size: Number of assembly accessions to batch submit to 'datasets'.

    Returns:
        Dictionary of core name and its associated assembly report
    """

    master_accn_list = list(assembly_accessions.values())
    combined_asm_reports = {}

    # Setting the number of combined accessions to query in a single call to datasets
    list_split = list(range(0, len(master_accn_list), batch_size))
    accn_subsample = [master_accn_list[ind : ind + batch_size] for ind in list_split]

    for accessions in accn_subsample:
        datasets_command = ["datasets", "summary", "genome", "accession"] + accessions

        # Make call to singularity datasets providing a multi-accession query:
        client_return = Client.execute(
            image=sif_image, command=datasets_command, return_result=True, quiet=True
        )

        raw_result = client_return["message"]

        ## Test what result we have obtained following execution of sif image and accession value
        # Returned a str, i.e. no datasets result obtained exited with fatal error
        if isinstance(raw_result, list):
            result = raw_result[0]
        else:
            result = raw_result

        if not isinstance(result, str):
            raise ValueError("Result obtained from datasets is not the expected format 'string'")
        if re.search("^FATAL", result):
            raise RuntimeError(f"Singularity image execution failed! -> '{result.strip()}'")
        # Returned a list, i.e. datasets returned a result to client.execute

        tmp_asm_dict = json.loads(result)
        if tmp_asm_dict["total_count"] >= 1:
            logging.info(f"Assembly report obtained for accession(s) {accessions}")

            batch_reports_json = tmp_asm_dict["reports"]
            for assembly_report in batch_reports_json:
                accession = assembly_report["accession"]
                asm_json_outfile = Path(download_directory, f"{accession}.asm_report.json")
                print_json(asm_json_outfile, assembly_report)
                # Save assembly report into master core<>report dict
                for core, accession_core in assembly_accessions.items():
                    if accession == accession_core:
                        combined_asm_reports[core] = assembly_report
        else:
            logging.warning(f"No assembly report found for accession(s) {accessions}. Exiting !")

    return combined_asm_reports


def extract_assembly_metadata(assembly_reports: Dict[str, dict]) -> Dict[str, ReportStructure]:
    """Function to parse assembly reports and extract specific key information on
    status and related fields.

    Args:
        assembly_reports: Key value pair of core_name : assembly report.

    Returns:
        Parsed assembly report meta (core, meta).
    """
    parsed_meta = {}

    for core, asm_report in assembly_reports.items():
        asm_meta_info = ReportStructure()

        # Mandatory meta key parsing:
        asm_meta_info["Asm accession"] = asm_report["accession"]
        asm_meta_info["Asm name"] = asm_report["assembly_info"]["assembly_name"]
        asm_meta_info["Assembly type"] = asm_report["assembly_info"]["assembly_type"]
        asm_meta_info["Asm status"] = asm_report["assembly_info"]["assembly_status"]
        asm_meta_info["Species Name"] = asm_report["organism"]["organism_name"]
        asm_meta_info["Taxon ID"] = asm_report["organism"]["tax_id"]

        ## Non-mandatory meta key parsing:
        assembly_meta_keys = asm_report["assembly_info"].keys()
        organism_keys = asm_report["organism"].keys()

        # check for genome_notes:
        if "genome_notes" in assembly_meta_keys:
            complete_notes = ", ".join(asm_report["assembly_info"]["genome_notes"])
            asm_meta_info["Asm notes"] = complete_notes
        else:
            asm_meta_info["Asm notes"] = "NA"

        # check for biosample:
        if "biosample" in assembly_meta_keys:
            asm_meta_info["Asm last updated"] = asm_report["assembly_info"]["biosample"]["last_updated"]
        else:
            asm_meta_info["Asm last updated"] = "NA"

        # check for paired assembly:
        if "paired_assembly" in assembly_meta_keys:
            asm_meta_info["Paired assembly"] = asm_report["assembly_info"]["paired_assembly"]["accession"]
        else:
            asm_meta_info["Paired assembly"] = "NA"

        # check for isolate/strain type:
        if "infraspecific_names" in organism_keys:
            organism_type_keys = asm_report["organism"]["infraspecific_names"].keys()
            if "isolate" in organism_type_keys:
                asm_meta_info["Isolate/Strain"] = asm_report["organism"]["infraspecific_names"]["isolate"]
            elif "strain" in organism_type_keys:
                asm_meta_info["Isolate/Strain"] = asm_report["organism"]["infraspecific_names"]["strain"]
            else:
                asm_meta_info["Isolate/Strain"] = "NA"
        else:
            asm_meta_info["Isolate/Strain"] = "NA"

        parsed_meta[core] = asm_meta_info

    return parsed_meta


def generate_report_tsv(
    parsed_asm_reports: Dict,
    outfile_prefix: str,
    query_type: str,
    output_directory: PathLike = Path(),
) -> None:
    """Generate and write the assembly report to a TSV file.

    Args:
        parsed_asm_reports: Parsed assembly report meta
        output_directory: Path to directory where output TSV is stored.
        query_type: Type of query core_db or accession
        output_directory: Directory to store report TSV
    """

    if not parsed_asm_reports:
        return

    tsv_outfile = Path(output_directory, f"{outfile_prefix}.tsv")

    header_list = list(ReportStructure().keys())
    header_list = [query_type] + header_list

    with open(tsv_outfile, "w+") as tsv_out:
        writer = csv.writer(tsv_out, delimiter="\t", lineterminator="\n")
        writer.writerow(header_list)

        for core, report_meta in parsed_asm_reports.items():
            final_asm_report = [core] + list(report_meta.values())
            writer.writerow(final_asm_report)
        tsv_out.close()


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input_cores",
        type=Path,
        required=False,
        default=None,
        help="List of ensembl core database(s) names to retrieve query accessions",
    )
    input_group.add_argument(
        "--input_accessions",
        type=Path,
        required=False,
        default=None,
        help="List of assembly INSDC query accessions",
    )
    parser.add_argument_dst_path(
        "--download_dir",
        default="Assembly_report_jsons",
        help="Folder where the assembly report JSON file(s) are stored",
    )
    parser.add_argument_dst_path(
        "--assembly_report_prefix",
        default="AssemblyStatusReport",
        help="Prefix used in assembly report TSV output file.",
    )
    parser.add_argument(
        "--host",
        type=str,
        required=False,
        help="Server hostname (fmt: mysql-ens-XXXXX-YY); required with '--input_cores'",
    )
    parser.add_argument(
        "--port", type=str, required=False, help="Server port (fmt: 1234); required with '--input_cores'"
    )
    parser.add_argument(
        "--datasets_version_url",
        type=str,
        required=False,
        metavar="URL",
        help="datasets version: E.g. docker://ensemblorg/datasets-cli:latest",
    )
    parser.add_argument(
        "--cache_dir",
        type=Path,
        required=False,
        default="$NXF_SINGULARITY_CACHEDIR",
        metavar="SINGULARITY_CACHE",
        help="Custom path to user generated singularity container housing ncbi tool 'datasets'",
    )
    parser.add_argument(
        "--datasets_batch_size",
        type=int,
        required=False,
        default=100,
        metavar="BATCH_SIZE",
        help="Number of accessions requested in one query to datasets",
    )

    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()

    init_logging_with_args(args)

    # Set and create dir for download files
    args.download_dir.mkdir(parents=True, exist_ok=True)

    # Set input file and determine if proper parameterization options are defined.
    user_query_file = check_parameterization(args.input_cores, args.input_accessions, args.host, args.port)

    ## Parse and store cores/accessions from user input query file
    with user_query_file.open(mode="r") as f:
        query_list = f.read().splitlines()

    ## Parse singularity setting and define the SIF image for 'datasets'
    datasets_image = singularity_image_setter(args.cache_dir, args.datasets_version_url)

    ## Get accessions on cores list or use user accession list directly
    connection_url = URL.create(
        "mysql",
        host=args.host,
        port=args.port,
        username="ensro",
    )
    query_accessions, query_type = resolve_query_type(
        query_list, connection_url, args.input_cores, args.input_accessions
    )

    # Datasets query implementation for one or more batched accessions
    assembly_reports = datasets_asm_reports(
        datasets_image, query_accessions, args.download_dir, args.datasets_batch_size
    )

    # Extract the key assembly report meta information for reporting status
    key_assembly_report_meta = extract_assembly_metadata(assembly_reports)

    # Produce the finalized assembly status report TSV from set of parsed 'datasets summary report'
    generate_report_tsv(key_assembly_report_meta, args.assembly_report_prefix, query_type, args.download_dir)
