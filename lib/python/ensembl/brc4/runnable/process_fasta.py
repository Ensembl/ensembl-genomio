#!env python3

import os
import eHive
import gzip
import shutil

class process_fasta(eHive.BaseRunnable):

    def run(self):
        genome_data = self.param('genome_data')
        file_name = self.param("file_name")
        file_path = self.param(file_name)
        work_dir = self.param('work_dir')

        if not self.param_exists(file_name):
            return

        accession = genome_data["assembly"]["accession"]

        # Set and create dedicated dir for download
        work_dir = os.path.join(work_dir, accession)
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final file name
        new_file_name = file_name + ".fa"

        # Final path
        head, tail = os.path.split(file_path)
        final_path = os.path.join(work_dir, new_file_name)

        # Uncompress or copy
        if file_path.endswith(".gz"):
            with open(final_path, "wb") as out_fasta:
                with gzip.open(file_path, "rb") as in_fasta:
                    out_fasta.write(in_fasta.read())
        else:
            shutil.rename(file_path, final_path)

        # No other operation
        self.dataflow({ file_name : final_path }, 2)

