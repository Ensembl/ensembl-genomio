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
"""Survey available genome assemblies from a taxonomic set of queries using NCBI's datasets CLI tool."""

from enum import Enum
import json
import logging
import os
import pandas as pd
from pathlib import Path
import re

from spython.main import Client


import ensembl.io.genomio
from ensembl.io.genomio.assembly.status import singularity_image_setter
# from ensembl.io.genomio.assembly.status import DATASETS_SINGULARITY
from ensembl.io.genomio.assembly.status import fetch_datasets_reports
from ensembl.io.genomio.utils.json_utils import print_json

from ensembl.ncbi_datasets.taxonomy_report import TaxonomyReports
from ensembl.ncbi_datasets.genome_summary import GenomeAssemblyReports

from ensembl.utils import StrPath
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args

# class DatasetsPackage(Enum):
#         GENE = "gene"
#         GENOME = "genome"
#         TAXONOMY = "taxonomy"
#         VIRUS = "virus"

def fetch_datasets_data_package(sif_image: Client, data_package: str, queries: str, download_dir: StrPath) -> dict[str, TaxonomyReports]:
    """A generic datasets-cli caller, to fetch data packages from NCBI and return the raw combined result.

    Args:
        sif_image: Instance of `Client.loaded()` singularity image.
        data_package: 
        queries:
        download_dir: Directory path to store parsed datasets-cli reports
    """

    # # Type checking
    # print(type(data_package))
    # if not isinstance(data_package, DatasetsPackage):
    #     raise TypeError(f"data package '{data_package}' must be an instance of 'DatasetsPackage' Enum")
    
    # base_datasets_cmd = ["datasets", "summary", f"{data_package}"]
    datasets_command = ["datasets", "summary", "taxonomy", "taxon", f"{queries}"]

    # Setting the number of combined accessions to query in a single call to datasets (50 max queries)
    sample_limit = 100
    list_split = list(range(0, len(queries), sample_limit))
    query_subsamples = [list_split[ind : ind + sample_limit] for ind in list_split]
    
    print(f"command: {datasets_command}")
    for sub_sample in query_subsamples:
        # Make call to singularity datasets providing a multi-accession query
        client_return = Client.execute(
            image=sif_image, command=datasets_command, return_result=True, quiet=True
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
            logging.warning(f"No datasets report found for query: {sub_sample}")
            continue

        logging.info(f"Taxonomy report obtained for queries(s) {queries}")

        taxonomy_reports: dict = {}
        batch_reports_json = tmp_asm_dict["reports"]
        for report in batch_reports_json:
            taxonReport = TaxonomyReports(json.dumps(report))
            species_name = taxonReport.species_name.replace(" ", "_")
            taxon_id = taxonReport.ncbi_taxon_id
            sci_taxon_name = f"{species_name}-{taxon_id}"
            taxonomy_reports[f"{sci_taxon_name}"] = taxonReport
            asm_json_outfile = Path(download_dir, f"{sci_taxon_name}.taxon-report.json")
            print_json(asm_json_outfile, report)

    return taxonomy_reports


            # asm_json_outfile = Path(download_directory, f"{accession}.json")
            # print_json(asm_json_outfile, assembly_report)
            # # Save assembly report into source key<>report dict
            # for src_key, accession_core in assembly_accessions.items():
            #     if accession == accession_core:
            #         combined_asm_reports[src_key] = assembly_report

# def fetch_taxonomy_reports(download_directory: StrPath):
#     """
#     TODO:
    
#     Args:
#         download_directory: Directory path to store assembly report JSON files.
#     """


def survey_genome_setup_from_tsv (src_file: StrPath) -> str:
    """Function to parse input TSV and partition appropriate survey information required for generating reports.
    
    Args:
        src_file: Path to file with one line per query taxonomic species/groups
    
    Returns:
        queries: list of taxon_ids for which to search
    """
    queries:list = []
    taxonomy_df = pd.read_csv(f"{src_file}", sep='\t')
    # Iterate over pandas dataframe
    for data_row in taxonomy_df.itertuples(index=False):
        taxon_id = data_row[0]
        queries.append(str(taxon_id))
    print(queries)
    taxon_list = ",".join(queries)
    
    return taxon_list
        # name = data_row[1]
        # accession = data_row[2]
        # contact = data_row[3]
        # print(f"working on :taxon ID {taxon_id}")
        # print(f"working on :name {name}")
        # print(f"working on :accession {accession}")
        # print(f"working on :contact {contact}")

    

def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_argument("--input", required=True, help="Taxonomic input TSV used to construct queries")
    parser.add_argument_dst_path(
        "--reports_dir",
        default=Path("assembly_report_jsons"),
        help="path to folder where the assembly report JSON files is stored",
    )
    parser.add_argument_dst_path(
        "--download_dir", default=Path.cwd(), help="Folder where the data will be downloaded"
    )
    parser.add_argument(
        "--datasets_version_url",
        type=str,
        metavar="URL",
        default="docker://ensemblorg/datasets-cli:latest",
        help="datasets version, e.g. docker://ensemblorg/datasets-cli:latest",
    )
    parser.add_argument_src_path(
        "--cache_dir",
        default=Path(os.environ.get("NXF_SINGULARITY_CACHEDIR", "")),
        metavar="SINGULARITY_CACHE",
        help="folder path to user generated singularity container housing NCBI tool 'datasets'",
    )
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    # Parse input TSV for taxon queries
    taxon_queries = survey_genome_setup_from_tsv(args.input)
    
    # Parse singularity setting and define the SIF image for 'datasets'
    datasets_image = singularity_image_setter(args.cache_dir, args.datasets_version_url)

    # fetch_datasets_data_package(datasets_image, "taxonomy", taxon_queries, args.download_dir)
    fetch_datasets_data_package(datasets_image, "taxonomy", taxon_queries, args.download_dir)