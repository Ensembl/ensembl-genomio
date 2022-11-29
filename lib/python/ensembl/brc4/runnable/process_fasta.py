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



from functools import partial
import gzip
from mimetypes import guess_type
from pathlib import Path

from Bio import SeqIO
import eHive


class process_fasta(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                "exclude_seq_regions": [],
                "peptide" : False,
                }

    def run(self):
        file_name = self.param("file_name")
        file_path = self.param(file_name)
        work_dir = Path(self.param('work_dir'))
        seqr_to_exclude = self.param("exclude_seq_regions")
        
        if self.param("peptide"):
            genbank_path = self.param('in_genbank')
            to_exclude = self.peptides_to_exclude(genbank_path, seqr_to_exclude)
        else:
            to_exclude = seqr_to_exclude

        if not self.param_exists(file_name):
            return

        if not work_dir.is_dir():
            work_dir.mkdir(parents=True)

        # Final file name
        new_file_name = file_name + ".fa"

        # Final path
        final_path = work_dir / new_file_name

        # Copy and filter
        records = []
        encoding = guess_type(file_path)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        
        with _open(file_path) as in_fasta:
            for record in SeqIO.parse(in_fasta, "fasta"):
                if record.id in to_exclude:
                    print("Skip record %s" % record.id)
                else:
                    records.append(record)

        with final_path.open("w") as out_fasta:
            SeqIO.write(records, out_fasta, "fasta")

        # No other operation
        self.dataflow({ file_name : str(final_path) }, 2)
    
    def peptides_to_exclude(self, genbank_path, seqr_to_exclude) -> dict():
        """
        Extract peptide IDs from a genbank file that are in a given list of seq regions
        """
        encoding = guess_type(genbank_path)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        
        peptides_to_exclude = dict()
        with _open(genbank_path) as in_genbank:
            for record in SeqIO.parse(in_genbank, "genbank"):
                if record.id in seqr_to_exclude:
                    print("Skip sequence %s" % record.id)
                    for feat in record.features:
                        if feat.type == "CDS":
                            if "protein_id" in feat.qualifiers:
                                feat_id = feat.qualifiers["protein_id"]
                                peptides_to_exclude[feat_id[0]] = True
                            else:
                                raise Exception("Peptide without peptide ID %s" % feat)
        
        return peptides_to_exclude
