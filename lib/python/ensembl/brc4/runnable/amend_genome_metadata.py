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


import os
from pathlib import Path
import re
from typing import Any, Dict, List, Tuple
import eHive
import gzip
import csv

from Bio import SeqIO
from ensembl.brc4.runnable.utils import Utils


class MissingDataError(Exception):
    """Used if some data is missing from the report file."""
    def __init__(self, report_path: str, accession: str, msg):
        report_msg = f"Can't get data for {accession} in report {report_path}"
        if msg:
            report_msg = f"{report_msg}: {msg}"
        self.msg = report_msg


class amend_genome_metadata(eHive.BaseRunnable):
    """Runnable to add things to the genome metadata, if any.

    Params:
        genome_data: Genome data following the schema/genome_schema.json.
        work_dir: Directory where the new file will be created.
        report: Path to the INSDC/RefSeq sequences report to parse.
        gbff: Path to the INSDC/RefSeq gbff file to parse.
        brc4_mode: Activate BRC4 mode (default).
    
    Dataflows:
        2: One dict with 2 keys:
            metadata_type: 'genome'.
            metadata_json: Path to the new genome json.
    
    Notes:
        BRC4 mode does the following:
            NOTHING FOR NOW
    """

    def param_defaults(self):
        return {
            "brc4_mode": True,
        }

    def run(self):
        genome_data = self.param('genome_data')
        work_dir = self.param('work_dir')
        report_path = self.param('report')
        gbff_path = self.param('gbff')
        brc4_mode = self.param('brc4_mode')

        # Create dedicated work dir
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final file name
        metadata_type = "genome"
        new_file_name = f"{metadata_type}_amended.json"
        final_path = Path(work_dir, new_file_name)
        
        use_refseq = self.param('accession').startswith("GCF_")
        
        # Get additional sequences in the assembly but not in the data
        additions = self.get_additions(report_path, gbff_path)

        if additions:
            genome_data["added_seq"] = {"region_name": additions}

        # Print out the file
        Utils.print_json(final_path, genome_data)

        # Flow out the file and type
        output = {
            "metadata_type": metadata_type,
            "metadata_json": str(final_path)
        }
        self.dataflow(output, 2)
    
    def get_additions(self, report_path: str, gbff_path: str) -> List[str]:
        """Returns all `seq_regions` that are mentioned in the report but that are not in the data.
        
        Args:
            report_path: Path to the report file.
            gbff_path: Path to the GBFF file.
            
        """
        gbff_regions = set(self.get_gbff_regions(gbff_path))
        report_regions = set(self.get_report_regions(report_path))

        additions = list(report_regions.difference(gbff_regions))
        return additions

    def get_gbff_regions(self, gbff_path: str) -> List[str]:
        """Returns the `seq_region` data from the GBFF file.

        Args:
            gbff_path: Gbff file path to use.

        """
        if not gbff_path:
            return []
        
        seq_regions = []
        _open = gbff_path.endswith(".gz") and gzip.open or open
        with _open(gbff_path, 'rt') as gbff_file:

            for record in SeqIO.parse(gbff_file, "genbank"):
                seq_regions.append(record.id)
                    
        return seq_regions
    
    def report_to_csv(self, report_path: str) -> Tuple[str, dict]:
        """Returns an assembly report as a CSV string, and the head metadata as a dict.

        Args:
            report_path: Path to a `seq_region` file from INSDC/RefSeq.

        """

        _open = report_path.endswith(".gz") and gzip.open or open
        with _open(report_path, 'rt') as report:
            data = ""
            metadata = {}
            last_head = ""
            for line in report:
                # Ignore header
                if line.startswith("#"):
                    # Get metadata values if possible
                    match = re.search("# (.+?): (.+?)$", line)
                    if match:
                        metadata[match.group(1)] = match.group(2)
                    last_head = line
                    continue
                else:
                    if last_head:
                        data += last_head[2:].strip() + "\n"
                        last_head = ""
                    data += line
            return data, metadata
    
    def get_report_regions(self, report_path: str) -> List[str]:
        """Returns a list of `seq_region` names from the report file.

        Args:
            report_path: Path to the seq_regions report from INSDC/RefSeq.

        """

        # Get the report in a CSV format, easier to manipulate
        report_csv, _ = self.report_to_csv(report_path)
        
        # Feed the csv string to the CSV reader
        reader = csv.DictReader(report_csv.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE)
        
        # Create the seq_regions
        seq_regions = []
        for row in reader:
            refseq_name = row["RefSeq-Accn"]
            genbank_name = row["GenBank-Accn"]
            name = ''
            if genbank_name and genbank_name != "na":
                name = genbank_name
            elif refseq_name and refseq_name != "na":
                name = refseq_name
            if name:
                seq_regions.append(name)
        
        return seq_regions
