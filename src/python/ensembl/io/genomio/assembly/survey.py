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
# from collections import defaultdict
import json
import logging
import os
import pandas as pd
from pathlib import Path
import re
import sys

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


class DatasetsPackage(Enum):
    GENE = "gene"
    GENOME = "genome"
    TAXONOMY = "taxonomy"
    VIRUS = "virus"


def fetch_datasets_data_package(
    sif_image: Client,
    data_package: DatasetsPackage,
    queries: list,
    download_dir: StrPath,
    extra_params: str | None,
) -> dict[str, TaxonomyReports]:
    """A generic datasets-cli caller, to fetch data packages from NCBI and return the raw combined result.

    Args:
        sif_image: Instance of `Client.loaded()` singularity image.
        data_package: Instance of DatasetsPackage(Enum): [taxonomy, genome, gene, virus]
        queries: User specific queries (accessions, taxon_ids, gene_ids etc.)
        download_dir: Directory path to store parsed datasets-cli reports
        extra_params: Additional params to append to dataset-cli command call (flags, positional arguments etc).
    """

    # # Type checking
    if not isinstance(data_package, Enum):
        raise TypeError(f"data package '{data_package}' must be an instance of 'DatasetsPackage' Enum")

    if data_package.value == "taxonomy":
        base_command = ["datasets", "summary", f"{data_package.value}", "taxon"]
        ReportClass = TaxonomyReports
    elif data_package.value == "genome":
        base_command = ["datasets", "summary", f"{data_package.value}", "accession"]
        ReportClass = GenomeAssemblyReports
    # elif data_package.value == "gene": TODO
    # elif data_package.value == "virus": TODO

    # Setting the number of combined accessions to query in a single call to datasets (50 max queries)
    batch_size_limit = 20  # CHANGE VALUE AFTER DEV!!
    list_split = list(range(0, len(queries), batch_size_limit))
    query_subsamples = [queries[ind : ind + batch_size_limit] for ind in list_split]

    for sub_sample in query_subsamples:
        # Make call to singularity datasets providing a multi-accession query
        datasets_command = base_command + sub_sample
        if extra_params:
            datasets_command += extra_params.split(" ")
        # print(f"Submitting datasets-cli command: '{datasets_command}'")
        client_return = Client.execute(
            image=sif_image, command=datasets_command, return_result=True, quiet=True
        )
        raw_result = client_return["message"]

        # Test what result we have obtained following execution of sif image and accession value
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

        # ingest report(s) into JSON dict
        json_slurped_reports = json.loads(result)
        if not json_slurped_reports["total_count"]:
            logging.warning(f"No datasets report(s) found for query: {sub_sample}")
            continue

        logging.info(f"Datasets report (type: {type(ReportClass)} obtained for queries(s) {sub_sample}")

        # Parse reports to dump to file and store parsed_report_object
        parsed_ncbi_reports: dict = {}
        batch_reports_json = json_slurped_reports["reports"]
        for report in batch_reports_json:
            json_report = ReportClass(json.dumps(report))
            # Formatted data package specific report|file name
            report_name, file_suffix = prepare_file_output(data_package.name, json_report)
            parsed_ncbi_reports[f"{report_name}"] = report
            asm_json_outfile = Path(download_dir, f"{report_name}-{file_suffix}")
            print_json(asm_json_outfile, report)
    return parsed_ncbi_reports


def prepare_file_output(datapackage_type: str, report: dict) -> str:
    """
    Args:
        datapackage:

    """
    # print(f"report has type: {type(report)}")

    if datapackage_type == "TAXONOMY":
        species_name = report.species_name.replace(" ", "_")
        taxon_id = report.ncbi_taxon_id
        file_prefix = f"{species_name}"
        file_suffix = "taxonReport.json"
    if datapackage_type == "GENOME":
        accession = report.accession
        taxon_id = report.organism["taxId"]
        file_prefix = f"{accession}"
        file_suffix = "genomeReport.json"
    else:
        print("REPORT FILE NAME FORMATER EXITED")
        sys.exit()

    fmt_outfile_name = f"{file_prefix}-{taxon_id}"
    # print(f"PREFIX = {file_prefix}, SUFFIX = {file_suffix}, FILENAME = {fmt_outfile_name}")
    return fmt_outfile_name, file_suffix


# def survey_setup_from_tsv(src_file: StrPath) -> dict:
def parse_query_tsv(src_file: StrPath) -> dict:
    """Function to parse input TSV and partition appropriate survey information required for generating reports.

    Args:
        src_file: Path to file with one line per query taxonomic|genome species/groups

    Returns:
        Dictionary of query lists of accessions and/or taxon IDs and/or species names).
    """
    parsed_queries: dict[list] = {"Accessions": [], "Taxon_ids": [], "Species": []}
    query_count = 0
    parsed_accession, parsed_taxonid, parsed_species = [], [], []

    # Iterate over pandas dataframe
    taxonomy_df = pd.read_csv(f"{src_file}", sep="\t")
    taxonomy_df.fillna("NA", inplace=True)
    for data_row in taxonomy_df.itertuples(index=False):
        if data_row[0] != "NA":
            taxon_id = str(int(data_row[0]))
        else:
            taxon_id = data_row[0]

        # Check taxonID and accessions present and not NA:
        taxon_match = re.match(r"^[0-9]+", taxon_id)

        ## ALSO NEED TO ADD REGEX FOR BIOPROJECT ACCESSION TODO
        accession_match = re.match(r"^(GC[AF])_([0-9]{3})([0-9]{3})([0-9]{3})(\.[0-9]+)?$", data_row[2])
        sp_name = data_row[1]
        opt_params = data_row[3]
        contact = data_row[4]

        # Best case, use accession if present
        if accession_match:
            parsed_accession.append(
                {"insdc_accession": accession_match.group(0), "opt_params": opt_params, "contact": contact}
            )
        # TaxonID present, accession missing:
        elif taxon_match and not accession_match:
            parsed_taxonid.append(
                {"taxon_id": taxon_match.group(0), "opt_params": opt_params, "contact": contact}
            )
        # Both taxonID and accession is missing:
        elif sp_name != "NA":
            parsed_species.append({"organism": sp_name, "opt_params": opt_params, "contact": contact})
        else:
            logging.warning(
                f"Unable to fully parse row {str(data_row)} in TSV, no useful species, genome or taxon ID info present"
            )
        query_count += 1
    logging.info(f"Input TSV parsed, parsed n={query_count} queries for further processing")
    logging.info(
        f"Breakdown of queries: Accession based n=({len(parsed_accession)}),TaxonID based n=({len(parsed_taxonid)}) and Species based n=({len(parsed_species)})"
    )
    
    parsed_queries["Accessions"] = parsed_accession
    parsed_queries["Taxon_ids"] = parsed_taxonid
    parsed_queries["Species"] = parsed_species
    # print(parsed_queries)
    return parsed_queries

def setup_genome_survey(sif_cache_dir: Path | None, datasets_version: str | None, input_queries: StrPath, out_dir: StrPath | None) -> None:
    """Main function to download and process genome reports downloaded using datasets-cli
    
    Args:
        cache_dir:
        datasets_version: 
        input_queries:
        download_dir:
    
    Returns:
        None.... for now !
    """

    # Parse input TSV for taxon queries
    parsed_queries: dict = {}
    parsed_queries = parse_query_tsv(input_queries)

    # Parse singularity setting and define the SIF image for 'datasets'
    datasets_image: Client = singularity_image_setter(sif_cache_dir, datasets_version)

    # Iterate over parsed queries and generate raw report JSONs
    for query_type, queries in parsed_queries.items():

        # Download exact genome report from specific accession
        if query_type == "Accessions":
            # print(f"Iterating over accession based queries -> {queries}")
            # Group accession query dict(list) by optional params in order to batch submit accessions accurately to datasets-cli
            opt_params_group: list = (pd.DataFrame(queries)
                        .groupby(['opt_params'])
                        .agg(list)
                        .reset_index()
                        .to_dict('records'))
            # Parse optional group dicts and submit each set of accession queries with the same input optional params
            for grouped_params in opt_params_group:
                accessions=grouped_params["insdc_accession"]
                datasets_params=grouped_params["opt_params"]
                # print(f"We have the following accessions: {accessions} to be searched via this param: {datasets_param}")
                #Call datasets cli with optional params
                if datasets_params != "NA":
                    logging.info(f"Calling datasets on accession based queries {accessions} + with user specific datasets-cli parameters {datasets_params}.")
                    fetch_datasets_data_package(datasets_image, DatasetsPackage.GENOME, accessions, out_dir, datasets_params)
                #call datasets cli with None for optional params
                else:
                    logging.info(f"Calling datasets on accession based queries {accessions}. No additional datasets-cli params specified.")
                    fetch_datasets_data_package(datasets_image, DatasetsPackage.GENOME, accessions, out_dir, None)
        # # Download set of genome reports from taxon ID
        # elif query_type == "Taxon_ids":
        #     print(f"Iterating over taxon ID based queries -> {queries}")
        # # Download set of genome reports from species name
        # elif query_type == "Species":
        #     print(f"Iterating over species based queries -> {queries}")

def main() -> None:
    """Module's entry-point."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_argument("--input", required=True, help="Input TSV used to construct queries")
    parser.add_argument_dst_path(
        "--output_dir",
        default=Path("assembly_report_jsons"),
        help="Path to folder where the assembly report JSON files is stored",
    )
    parser.add_argument_dst_path(
        "--download_dir", default=Path.cwd(), help="Path to store downloaded genome report JSONs"
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

    setup_genome_survey(args.cache_dir, args.datasets_version_url, args.input, args.download_dir)
