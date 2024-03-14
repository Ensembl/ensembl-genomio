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
"""Record the assembly status for a set of INSDC accessions using ncbi 'datasets' tool"""

__all__ = [
    "fetch_asm_accn",
    "datasets_asm_report",
]

import json
import os
from os import PathLike, getcwd
from pathlib import Path
import re
import logging

from spython.main import Client

from ensembl.io.genomio.utils.json_utils import print_json
from ensembl.database import DBConnection as dbc
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

DATASETS_SINGULARITY = {
    "datasets_version_url": "library://lcampbell/ensembl-genomio/ncbi-datasets-v16.6.0:latest",
}

# class ImageDoesNotExist(Exception):
#     """When spython client.pull fails to locate a specified sif image"""

# def fetch_datasets_singularity(image_url: str, image_dl_path: str) -> Client.instance:
#     """Obtain the requested image container from singularity hub using spython package

#     Args:
#         image_url:
#         image_dl_path:

#     Returns:
#         container_image: Client.instance of singularity sif image loaded
#     """

#     # Pull container image to specified download path, if exists do nothing
#     datasets_image = Client.pull(image_url,
#                     stream=False,
#                     pull_folder=image_dl_path,
#                     quiet=True)

#     return datasets_image
# except ImageDoesNotExist:
#     logging.warning("While pulling library image: error fetching image: \
#                     image does not exist in the library")
# else:
#     print(type(datasets_image))
#     print(datasets_image)


def fetch_asm_accn(database_names: list, server_details: str) -> dict:
    """Obtain the associated INSDC accession [meta.assembly.accession] given a set of core(s) names
    and a MYSQL server host.

    Args:
        cores: Set of names for one or more core databases
        server_details: MYSQL host server name and port [mysql-ens-(NAME:Port)]

    Returns:
        Dict of core name(s) (key) and its INSDC assembly.accession (value)
    """

    core_accn_meta = {}
    core_list_count = len(database_names)
    count_accn_found = 0

    for core in database_names:
        db_connection_url = f"{server_details}{core}"
        db_connection = dbc(f"{db_connection_url}")
        qry_result = db_connection.execute(
            'SELECT meta_value FROM meta WHERE meta_key = "assembly.accession";'
        ).fetchall()

        if qry_result is None:
            logging.warning(f"We have no accession on core: {core}")
        elif len(qry_result) == 1:
            count_accn_found += 1
            asm_accession = qry_result.pop()[0]
            print(f"{core} -> assembly.accession[{asm_accession}]")
            core_accn_meta[core] = asm_accession
        else:
            logging.warning(f"Core {core} Has {len(qry_result)} assembly.accessions")

    logging.info(f"From initial input cores ({core_list_count}), obtained ({count_accn_found}) accessions")

    return core_accn_meta


# def single_datasets_asm_report(sif_image: str, assembly_accessions: dict, download_directory: PathLike) -> None:
#     """Obtain one single assembly report JSON per assembly accession,
#     i.e. make individual since accn query to datasets tool.

#     Args:
#         sif_image: Instance of Client.loaded singularity image.
#         assembly_accessions: Dict of core accessions.
#         download_directory: Dir path to store assembly report JSON files.
#     """

#     ## Attempt to run container to pull assembly reports
#     for accession in assembly_accessions.values():
#         asm_json_outfile = f"{getcwd()}/{download_directory}/{accession}.asm_report.json"

#         # Actually make call to singularity:
#         client_return = Client.execute(image=sif_image,
#                                            command=['datasets','summary' ,'genome' ,'accession', accession],
#                                            return_result=True,
#                                            quiet=True)

#         result = client_return["message"]

#         ## Test what result we have returned following execution of sif image and accession value
#         # Returned a str, i.e. no datasets result obtained exited with fatal error
#         if isinstance(result, str) and re.search("^FATAL", result):
#             logging.critical(f"Singularity image execution failed! -> '{result.strip()}'")
#         # Returned a list, i.e. datasets returned a result to client.execute
#         elif isinstance(result, list):
#             tmp_asm_dict = json.loads(result.pop(0))

#             if tmp_asm_dict['total_count'] == 0:
#                 logging.warning(f"No assembly report found for accession {accession}")
#             elif tmp_asm_dict['reports']:
#                 logging.info(f"Asm report obtained for accession [{accession}]")
#                 final_asm_json = tmp_asm_dict['reports'][0]
#                 print_json(Path(asm_json_outfile), final_asm_json)
#             else:
#                 print("Something not right!")
#         else:
#             logging.warning(f"Something not right while running datasets with singularity client {client_return}")


def datasets_asm_report(sif_image: str, assembly_accessions: dict, download_directory: PathLike, batch_size: int = 2) -> None:
    """Obtain multiple assembly report JSONs in one or more querys to datasets,
    i.e. make individual since accn query to datasets tool.

    Args:
        sif_image: Instance of Client.loaded singularity image.
        assembly_accessions: Dict of core accessions.
        download_directory: Dir path to store assembly report JSON files.
    """

    ## Attempt to run container to pull assembly reports
    master_accn_list = list(assembly_accessions.values())

    # Setting the number of combined accessions to query in a single call to datasets
    list_split = [i for i in range(0, len(master_accn_list), batch_size)]  ## Note best to use >=10
    accn_subsample = [master_accn_list[ind : ind + batch_size] for ind in list_split]

    for accessions in accn_subsample:
        datasets_command = ["datasets", "summary", "genome", "accession"] + accessions

        # Make call to singularity datasets providing a multi accn query:
        client_return = Client.execute(
            image=sif_image, command=datasets_command, return_result=True, quiet=True
        )

        result = client_return["message"]

        ## Test what result we have returned following execution of sif image and accession value
        # Returned a str, i.e. no datasets result obtained exited with fatal error
        if isinstance(result, str) and re.search("^FATAL", result):
            logging.critical(f"Singularity image execution failed! -> '{result.strip()}'")
        # Returned a list, i.e. datasets returned a result to client.execute
        elif isinstance(result, list):
            tmp_asm_dict = json.loads(result.pop(0))

            if tmp_asm_dict["total_count"] == 0:
                logging.warning(f"No assembly report found for accession(s) {accessions}")
            elif tmp_asm_dict["reports"]:
                logging.info(f"Asm report obtained for accession(s) [{accessions}]")
                batch_reports_json = tmp_asm_dict["reports"]
                for assembly_report in batch_reports_json:
                    accession = assembly_report["accession"]
                    asm_json_outfile = f"{getcwd()}/{download_directory}/{accession}.asm_report.json"
                    print_json(Path(asm_json_outfile), assembly_report)
            else:
                print("Somethings not right!")
        else:
            logging.warning(
                f"Something not right while running datasets with singularity client.execute {client_return}"
            )


# def classify_assembly_status(core_accessions: dict) -> None:
#     """Main function to pare set of core list and call ncbi datasets"""


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description="test")
    parser.add_argument_src_path("--input_cores", required=True, help="List of ensembl core db names")
    parser.add_argument_dst_path(
        "--download_dir", default=Path.cwd(), help="Folder where the assembly report JSON file(s) are stored"
    )
    parser.add_argument("--host", type=str, required=True, help="Server hostname (fmt: mysql-ens-XXXXX-YY)")
    parser.add_argument("--port", type=str, required=True, help="Server port")
    parser.add_argument(
        "--datasets_version_url",
        type=str,
        required=False,
        metavar="URL",
        help="Custom datasets version. E.g. library://lcampbell/ensembl-genomio/ncbi-datasets-v16.6.0:latest",
    )
    parser.add_argument(
        "--cache_dir",
        type=Path,
        required=False,
        metavar="SINGULARITY_CACHE",
        help="Custom path to user generated singularity container housing ncbi tool 'datasets'",
    )

    parser.add_log_arguments()
    args = parser.parse_args()

    init_logging_with_args(args)

    # Set and create dir for download files
    if not args.download_dir.is_dir():
        args.download_dir.mkdir(parents=True)

    ## Main starts here:
    cores_list = []
    try:
        with args.input_cores.open(mode="r") as f:
            cores_list = f.read().splitlines()
    except IOError as err:
        logging.critical(f"Unable to read from database list inputfile '{args.input_cores}' due to {err}.")

    # Set singularity cache dir from user defined path or use environment
    if args.cache_dir and args.cache_dir.is_dir():
        image_dl_path = Path(args.cache_dir)
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
    if args.datasets_version_url is None:
        container_url = DATASETS_SINGULARITY["datasets_version_url"]
        logging.info(f"Using default 'ncbi datasets' version '{container_url}'")
    else:
        container_url = args.datasets_version_url
        logging.info(f"Using user defined 'ncbi datasets' version '{container_url}'")

    ## Get accessions on cores
    server_details = f"mysql://ensro@{args.host}:{args.port}/"
    core_db_accessions = fetch_asm_accn(cores_list, server_details)

    # Pull or load pre-existing 'datasets' singularity container image.
    # datasets_image = fetch_datasets_singularity(container_url, image_dl_path)
    datasets_image = Client.pull(container_url, stream=False, pull_folder=image_dl_path, quiet=True)

    # ## Attempt to run container to pull assembly reports
    # Test case implementing one ncbi datasets query per accession(s)
    # single_datasets_asm_report(datasets_image, core_db_accessions, args.download_dir)

    # Datasets query implementation for one or more bacthed accessions
    datasets_asm_report(datasets_image, core_db_accessions, args.download_dir)
