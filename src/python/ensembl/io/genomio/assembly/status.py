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
    "datasets_asm_reports",
    "extract_assembly_metadata",
    "generate_report_tsv",
]

import csv
import json
import os
from os import PathLike, getcwd
from pathlib import Path
import re
import logging
from typing import Dict

from spython.main import Client

from ensembl.io.genomio.utils.json_utils import print_json
from ensembl.database import DBConnection as dbc
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

DATASETS_SINGULARITY = {
    "datasets_version_url": "library://lcampbell/ensembl-genomio/ncbi-datasets-v16.6.0:latest",
}

class ReportStructure(dict):
    """Dict setter class of key report meta information"""
    def __init__(self):
        dict.__init__(self)
        self.update({
            'Species Name' : '',
            'Taxon ID' : '',
            'Strain' : '',
            'Isolate' : '',
            'Isolate/Strain' : '',
            'Asm name' : '',
            'Assembly type' : '',
            'Asm accession' : '',
            'Paired assembly' : '',
            'Asm last updated' : '',
            'Asm status' : '',
            'Asm notes' : '',
        })

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

def datasets_asm_reports(sif_image: str, assembly_accessions: dict,
                        download_directory: PathLike, batch_size: int = 20) -> dict:
    """Obtain multiple assembly report JSONs in one or more querys to datasets,
    i.e. make individual since accn query to datasets tool.

    Args:
        sif_image: Instance of Client.loaded singularity image.
        assembly_accessions: Dict of core accessions.
        download_directory: Dir path to store assembly report JSON files.
        batch_size: Number of assembly accessions to batch submit to 'datasets'.
    
    Returns:
        Dictionary of core name and its assoicated assembly report
    """

    master_accn_list = list(assembly_accessions.values())
    combined_asm_reports = {}

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

                    # Save assembly report into master core<>report dict
                    for core, accession_core in assembly_accessions.items():
                        if accession == accession_core:
                            combined_asm_reports[core] = assembly_report
            else:
                print("Something is not right!")
        else:
            logging.warning(
                f"Something not right while running datasets with singularity client.execute {client_return}"
            )
    return combined_asm_reports

def extract_assembly_metadata(assembly_reports: Dict[str,dict]) -> Dict[str,ReportStructure]:
    """"Function to parse assembly reports and extract specific key information on 
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
        # asm_meta_info["Asm last updated"] = asm_report["assembly_info"]["biosample"]["last_updated"]
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
                asm_meta_info["Isolate"] = asm_report["organism"]["infraspecific_names"]["isolate"]
                asm_meta_info.pop("Strain")
                asm_meta_info.pop("Isolate/Strain")
            elif "strain" in organism_type_keys:
                asm_meta_info["Strain"] = asm_report["organism"]["infraspecific_names"]["strain"]
                asm_meta_info.pop("Isolate")
                asm_meta_info.pop("Isolate/Strain")
            else:
                asm_meta_info["Isolate/Strain"] = 'NA'
                asm_meta_info.pop("Strain")
                asm_meta_info.pop("Isolate")
        else:
            # elif ("strain" not in organism_type_keys) and ("isolate" not in organism_type_keys):
            asm_meta_info["Isolate/Strain"] = 'NA'
            asm_meta_info.pop("Strain")
            asm_meta_info.pop("Isolate")

        parsed_meta[core] = asm_meta_info

    return parsed_meta

def generate_report_tsv(parsed_asm_reports: dict,  
                        outfile_prefix: str, 
                        output_directoy: PathLike = Path(getcwd())
                        ) -> None:
    """Generate and write the assembly report to a TSV file
    
    Args:
        parsed_asm_reports: Parsed assembly report meta
        output_directoy: Path to directory where output TSV is stored.
    """

    tsv_outfile = f"{output_directoy}/{outfile_prefix}.tsv"

    header_list = list(ReportStructure().keys())
    header_list.remove('Strain')
    header_list.remove('Isolate')
    header_list = ["Core DB"] + header_list

    with open(tsv_outfile, "w+") as tsv_out:

        writer = csv.writer(tsv_out, delimiter="\t", lineterminator="\n")
        writer.writerow(header_list)

        for core, report_meta in parsed_asm_reports.items():
            final_asm_report = [core] + list(report_meta.values())
            writer.writerow(final_asm_report)
        tsv_out.close()

# def isolate_pertinent_status(parsed_asm_reports: dict) -> None:

# def classify_assembly_status(core_accessions: dict) -> None:
#     """Main function to pare set of core list and call ncbi datasets"""


def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(
        description="Track the assembly status of a set of input core(s) using NCBI 'datasets'")
    parser.add_argument_src_path("--input_cores", required=True, help="List of ensembl core db names")
    parser.add_argument_dst_path(
        "--download_dir",
        default=Path.cwd(),
        help="Folder where the assembly report JSON file(s) are stored"
    )
    parser.add_argument_dst_path(
        "--assembly_report_prefix",
        default="AssemblyStatusReport",
        help="Prefix used in assembly report TSV output file."
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
        default="$NXF_SINGULARITY_CACHEDIR",
        metavar="SINGULARITY_CACHE",
        help="Custom path to user generated singularity container housing ncbi tool 'datasets'",
    )
    parser.add_argument(
        "--datasets_batch_size", 
        type=int,
        required=False,
        default=20,
        metavar="BATCH_SIZE",
        help="Number of accessions requested in one query to datasets"
    )

    parser.add_log_arguments(add_log_file=True)
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
    datasets_image = Client.pull(container_url, stream=False, pull_folder=image_dl_path, quiet=True)

    # Datasets query implementation for one or more bacthed accessions
    assembly_reports = datasets_asm_reports(datasets_image, core_db_accessions, args.download_dir)

    # Extract the key assembly report meta information for reporting status
    key_asmreport_meta = extract_assembly_metadata(assembly_reports)

    generate_report_tsv(key_asmreport_meta, args.assembly_report_prefix, args.download_dir)

    ## Parse reports and email key assembly report status (not current)
